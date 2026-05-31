#include "circuit.h"
#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <climits>
#include <random>

using namespace std;
using namespace Qcircuit;

// ================================================================
// CSR 빌더: interactionGraph → METIS CSR 포맷
//   perturb_ratio > 0 이면 edge weight에 ±perturb_ratio 비율 노이즈 추가
//   → 동일 그래프에서 다양한 METIS 파티션 탐색 가능
// ================================================================
struct MetisCSR {
    vector<idx_t> xadj;
    vector<idx_t> adjncy;
    vector<idx_t> adjwgt;
    vector<idx_t> vwgt;
    idx_t nvtxs;
};

static MetisCSR build_csr(Graph& g, double perturb_ratio, mt19937& rng)
{
    MetisCSR csr;
    csr.nvtxs = (idx_t)g.nodeset.size();
    csr.xadj.reserve(csr.nvtxs + 1);
    csr.xadj.push_back(0);

    uniform_real_distribution<double> noise(1.0 - perturb_ratio,
                                            1.0 + perturb_ratio);

    for (auto& [id, node] : g.nodeset)
    {
        for (int eid : node.edges)
        {
            auto& e = g.edgeset[eid];
            int neighbor = (e.getsourceid() == id) ? e.gettargetid()
                                                    : e.getsourceid();
            double w = e.getweight();
            if (perturb_ratio > 0.0)
                w *= noise(rng);

            idx_t iw = (idx_t)max(1LL, llround(w * 1000.0));
            csr.adjncy.push_back((idx_t)neighbor);
            csr.adjwgt.push_back(iw);
        }
        csr.xadj.push_back((idx_t)csr.adjncy.size());
    }
    csr.vwgt.assign(csr.nvtxs, 1);
    return csr;
}

// ================================================================
// 모듈성 지표 Q
//   Q = (파티션 내부 edge weight 합) / (전체 edge weight 합)
//   - 높을수록 inter-QPU gate 적음 (좋은 파티션)
//   - 빈 파티션: 강한 패널티로 사실상 제외
//   - 균형 패널티: 크기 편차 반영
// ================================================================
static double compute_modularity(
    Graph& g,
    const vector<idx_t>& part,
    int num_parts)
{
    double intra_w = 0.0;
    double total_w = 0.0;

    for (auto& [eid, e] : g.edgeset)
    {
        int s = e.getsourceid();
        int t = e.gettargetid();
        double w = e.getweight();
        total_w += w;
        if (s < (int)part.size() && t < (int)part.size() &&
            part[s] == part[t])
            intra_w += w;
    }

    if (total_w < 1e-12) return 0.0;

    double Q = intra_w / total_w;

    // 파티션 크기 분석
    vector<int> sz(num_parts, 0);
    for (auto v : part)
        if (v >= 0 && v < num_parts) sz[v]++;

    // 빈 파티션: 강한 패널티
    int empty_count = 0;
    for (int c : sz) if (c == 0) empty_count++;
    if (empty_count > 0)
        return Q - 1.0 * empty_count;

    // 불균형 패널티 계산
    double mean = (double)part.size() / num_parts;
    double var  = 0.0;
    for (int c : sz) var += (c - mean) * (c - mean);
    double stddev = sqrt(var / num_parts);
    double cv = stddev / (mean + 1e-9);  // 변동계수

    // 단계적 패널티:
    //  - 경미한 불균형 (cv≤0.3): 작은 패널티
    //  - 심한 불균형 (cv>0.3): 강한 패널티
    int max_sz = *max_element(sz.begin(), sz.end());
    double imbalance = (double)max_sz / (mean + 1e-9);  // 최대/평균 비율
    double penalty;
    if (imbalance <= 1.5) {
        penalty = 0.05 * cv;
    } else if (imbalance <= 2.0) {
        penalty = 0.15 * cv;   // 중간
    } else {
        penalty = 0.40 * cv;   // 심한 불균형 강하게 억제
    }

    return Q - penalty;
}

// ================================================================
// 빈 파티션 post-processing 재균형
//   빈 파티션이 있으면 가장 큰 파티션에서 노드를 하나씩 이동
// ================================================================
static void rebalance_partition(
    vector<idx_t>& part, int num_parts, int nvtxs)
{
    vector<vector<int>> members(num_parts);
    for (int i = 0; i < nvtxs; i++)
        if (part[i] >= 0 && part[i] < num_parts)
            members[part[i]].push_back(i);

    for (int p = 0; p < num_parts; p++)
    {
        if (!members[p].empty()) continue;
        // 가장 큰 파티션에서 이동
        int biggest = 0;
        for (int q = 1; q < num_parts; q++)
            if (members[q].size() > members[biggest].size()) biggest = q;
        if (members[biggest].size() <= 1) break;

        int node = members[biggest].back();
        members[biggest].pop_back();
        members[p].push_back(node);
        part[node] = (idx_t)p;
    }
}

// ================================================================
// 단일 METIS 실행
//   perturb_ratio: 가중치 교란 비율 (0이면 원본 그대로)
//   실패시 기본 옵션(nullptr)으로 재시도
// ================================================================
static bool run_metis_once(
    Graph& g,
    int num_parts,
    int seed,
    double perturb_ratio,
    mt19937& rng,
    vector<idx_t>& part_out,
    idx_t& objval_out)
{
    MetisCSR csr = build_csr(g, perturb_ratio, rng);
    idx_t nvtxs  = csr.nvtxs;
    idx_t ncon   = 1;
    idx_t nparts = (idx_t)num_parts;

    part_out.assign(nvtxs, 0);
    objval_out = 0;

    idx_t opts[METIS_NOPTIONS];
    METIS_SetDefaultOptions(opts);
    opts[METIS_OPTION_SEED]  = (idx_t)((long long)seed * 1000003LL % 999983 + 1);
    opts[METIS_OPTION_NITER] = 20;

    int status = METIS_PartGraphKway(
        &nvtxs, &ncon,
        csr.xadj.data(), csr.adjncy.data(),
        csr.vwgt.data(), nullptr, csr.adjwgt.data(),
        &nparts, nullptr, nullptr,
        opts, &objval_out, part_out.data()
    );

    // 실패시 기본 옵션으로 재시도
    if (status != METIS_OK) {
        MetisCSR csr2 = build_csr(g, 0.0, rng);   // no perturbation
        idx_t nv2 = csr2.nvtxs;
        part_out.assign(nv2, 0); objval_out = 0;
        status = METIS_PartGraphKway(
            &nv2, &ncon,
            csr2.xadj.data(), csr2.adjncy.data(),
            csr2.vwgt.data(), nullptr, csr2.adjwgt.data(),
            &nparts, nullptr, nullptr,
            nullptr, &objval_out, part_out.data()
        );
    }

    return (status == METIS_OK);
}

// ================================================================
// Adaptive Partitioning (메인 함수)
//
// 핵심 전략:
//   1. 첫 몇 회 (PURE_TRIALS): perturbation 없이 순수 METIS
//   2. 이후 (PERTURB_TRIALS): edge weight에 점진적 노이즈 추가
//      → 다양한 파티션 탐색 (perturbation이 커질수록 더 다른 결과)
//   3. 포화 감지: PATIENCE 연속 미개선 시 조기 종료
//   4. 재균형: 빈 파티션 post-processing 보정
//   5. 최고 Q 파티션 선택
// ================================================================
void Qcircuit::QMapper::run_metis_partition(Graph& interactiongraph, int num_parts)
{
    idx_t nvtxs = (idx_t)interactiongraph.nodeset.size();

    if (nvtxs <= (idx_t)num_parts) {
        partition_result.resize(nvtxs);
        for (idx_t i = 0; i < nvtxs; i++)
            partition_result[i] = i % num_parts;
        cout << "[ADAPTIVE-METIS] Trivial partition (nvtxs=" << nvtxs
             << " <= num_parts=" << num_parts << ")\n";
        return;
    }

    int num_edges = (int)interactiongraph.edgeset.size();

    // ── 회로 복잡도별 파라미터 설정 ────────────────────────────
    //   단순 회로 (≤30 edges): perturbation 불필요, 1회
    //   중간   (≤100 edges): perturbation 조금, 20회
    //   복잡   (>100 edges): perturbation 적극, 80회
    int  MAX_TRIALS, PATIENCE, PURE_TRIALS;
    double MAX_PERTURB;

    if (num_edges <= 30) {
        MAX_TRIALS  = 1;
        PATIENCE    = 1;
        PURE_TRIALS = 1;
        MAX_PERTURB = 0.0;
    } else if (num_edges <= 100) {
        MAX_TRIALS  = 40;
        PATIENCE    = 10;
        PURE_TRIALS = 2;    // 빠르게 perturbation 시작
        MAX_PERTURB = 0.40; // 더 큰 교란 허용
    } else {
        MAX_TRIALS  = 80;
        PATIENCE    = 15;
        PURE_TRIALS = 3;
        MAX_PERTURB = 0.45;
    }

    double IMPROVE_THRESH = 5e-4;

    // ── 반복 실행 ────────────────────────────────────────────────
    double        best_Q   = -1e18;
    idx_t         best_obj = (idx_t)INT_MAX;
    vector<idx_t> best_part;
    int no_improve = 0;
    int trial      = 0;
    int success    = 0;

    // RNG: trial마다 다른 노이즈
    mt19937 rng(42);

    cout << "\n[ADAPTIVE-METIS] start:"
         << " nvtxs=" << nvtxs
         << " edges=" << num_edges
         << " max_trials=" << MAX_TRIALS
         << " patience=" << PATIENCE
         << " max_perturb=" << MAX_PERTURB
         << "\n";

    while (trial < MAX_TRIALS)
    {
        // perturbation 비율: 처음 PURE_TRIALS는 0, 이후 점진적 증가
        double perturb = 0.0;
        if (trial >= PURE_TRIALS) {
            // 선형 증가: PURE_TRIALS ~ MAX_TRIALS 구간
            double progress = (double)(trial - PURE_TRIALS) /
                              max(1, MAX_TRIALS - PURE_TRIALS);
            perturb = MAX_PERTURB * progress;
        }

        // seed: 다양성 확보 (큰 소수 간격)
        int seed = trial * 137 + 42;

        vector<idx_t> part;
        idx_t objval;
        bool ok = run_metis_once(interactiongraph, num_parts, seed,
                                 perturb, rng, part, objval);

        if (!ok) { trial++; continue; }
        success++;

        // 빈 파티션 재균형
        rebalance_partition(part, num_parts, (int)nvtxs);

        double Q = compute_modularity(interactiongraph, part, num_parts);

        if (best_part.empty() || Q > best_Q + IMPROVE_THRESH)
        {
            double prev_Q = best_Q;
            best_Q    = Q;
            best_obj  = objval;
            best_part = part;
            no_improve = 0;

            cout << "  [trial " << setw(3) << trial
                 << "] Q=" << fixed << setprecision(5) << Q
                 << " cut=" << objval
                 << " perturb=" << setprecision(3) << perturb;
            if (prev_Q > -1e17)
                cout << " Δ=+" << setprecision(5) << (Q - prev_Q);
            cout << "  *** BEST ***\n";
        }
        else {
            no_improve++;
            if (trial % 10 == 9)
                cout << "  [trial " << setw(3) << trial
                     << "] Q=" << fixed << setprecision(5) << Q
                     << " best=" << best_Q
                     << " no_imp=" << no_improve << "/" << PATIENCE
                     << "\n";
        }

        if (no_improve >= PATIENCE) {
            cout << "  [ADAPTIVE-METIS] saturated at trial=" << trial
                 << " (no improvement for " << PATIENCE << " trials)\n";
            break;
        }
        trial++;
    }

    // ── fallback ─────────────────────────────────────────────────
    if (best_part.empty()) {
        cerr << "[ADAPTIVE-METIS] All trials failed, round-robin fallback\n";
        best_part.resize(nvtxs);
        for (idx_t i = 0; i < nvtxs; i++)
            best_part[i] = i % num_parts;
        best_Q = 0.0;
    }

    partition_result = std::move(best_part);

    // ── 결과 요약 ────────────────────────────────────────────────
    vector<int> sz(num_parts, 0);
    for (auto v : partition_result)
        if (v >= 0 && v < num_parts) sz[v]++;

    cout << "[ADAPTIVE-METIS] DONE:"
         << " best_Q=" << fixed << setprecision(5) << best_Q
         << " cut=" << best_obj
         << " trials=" << trial+1 << "/" << MAX_TRIALS
         << " success=" << success << "\n"
         << "  partition sizes: [";
    for (int i = 0; i < num_parts; i++)
        cout << sz[i] << (i+1<num_parts ? "," : "]\n");
}


// ================================================================
// analyze_cross_partition_edges
// ================================================================
void Qcircuit::QMapper::analyze_cross_partition_edges(
    Graph& g, const vector<idx_t>& part, int num_qpu)
{
    vector<vector<double>> cross(num_qpu, vector<double>(num_qpu, 0.0));

    for (const auto& [eid, e] : g.edgeset)
    {
        int s = e.getsourceid();
        int t = e.gettargetid();
        if (s < 0 || t < 0 ||
            s >= (int)part.size() || t >= (int)part.size()) continue;

        int ps = part[s], pt = part[t];
        double w = e.getweight();
        if (ps == pt) continue;
        cross[ps][pt] += w;
        cross[pt][ps] += w;
    }
    cross_table = cross;

    double total_cross = 0.0;
    for (int i = 0; i < num_qpu; i++)
        for (int j = i+1; j < num_qpu; j++)
            total_cross += cross[i][j];

    cout << "\n========== Cross-Subgraph Weighted Sums ==========\n"
         << "  Total cross-QPU weight: " << fixed << setprecision(3)
         << total_cross << "\n";

    vector<tuple<int,int,double>> cl;
    for (int i = 0; i < num_qpu; i++)
        for (int j = i+1; j < num_qpu; j++)
            cl.push_back({i, j, cross[i][j]});
    sort(cl.begin(), cl.end(),
         [](auto& a, auto& b){ return get<2>(a) > get<2>(b); });

    for (auto& [i, j, val] : cl)
        if (val > 1e-9)
            cout << "  QPU" << i << "<->QPU" << j
                 << ": " << setprecision(3) << val << "\n";
}


// ================================================================
// reconstruct_info
// ================================================================
void Qcircuit::QMapper::reconstruct_info(Graph& interactiongraph, int num_parts)
{
    Subsets.clear();
    Subsets.resize(num_parts);

    for (auto& [id, node] : interactiongraph.nodeset)
    {
        int p = partition_result[id];
        if (p >= 0 && p < num_parts)
            Subsets[p].nodeset[id] = node;
    }

    for (auto& [eid, e] : interactiongraph.edgeset)
    {
        int s = e.getsourceid(), t = e.gettargetid();
        int ps = partition_result[s], pt = partition_result[t];
        if (ps == pt && ps >= 0 && ps < num_parts)
        {
            Subsets[ps].edgeset[eid] = e;
            Subsets[ps].nodeset[s].edges.insert(eid);
            Subsets[ps].nodeset[t].edges.insert(eid);
        }
    }

    cout << "\n===== InteractionGraph Subsets Complete =====\n";
    inter_Sub_edge.clear();
    cout << "\n===== Building inter_Sub_edge (Cross-Sub Edges) =====\n";

    for (auto& [eid, e] : interactiongraph.edgeset)
    {
        int s = e.getsourceid(), t = e.gettargetid();
        int ps = partition_result[s], pt = partition_result[t];
        if (ps != pt)
            inter_Sub_edge.push_back({ ps, { pt, e } });
    }
}


// ================================================================
// matching
// ================================================================
void Qcircuit::QMapper::matching()
{
    cout << "\n===== MATCHING START =====\n";
    int num_sub = (int)Subsets.size();
    if (num_sub == 0) { cout << "No Subsets found.\n"; return; }

    matching_info.clear();
    matching_info.resize(num_sub, { Graph(), -1 });
    for (int i = 0; i < num_sub; i++)
    {
        matching_info[i] = { Subsets[i], i };
        cout << "[MATCH] Sub" << i << " -> QPU" << i << "\n";
    }
    cout << "\n===== MATCHING DONE =====\n";
}