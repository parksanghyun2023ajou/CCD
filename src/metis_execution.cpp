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

struct MetisCSR {
    vector<idx_t> xadj, adjncy, adjwgt, vwgt;
    idx_t nvtxs;
};

static MetisCSR build_csr(Graph& g, double perturb_ratio, mt19937& rng)
{
    MetisCSR csr;
    csr.nvtxs = (idx_t)g.nodeset.size();
    csr.xadj.push_back(0);
    uniform_real_distribution<double> noise(1.0 - perturb_ratio, 1.0 + perturb_ratio);
    for (auto& [id, node] : g.nodeset) {
        for (int eid : node.edges) {
            auto& e = g.edgeset[eid];
            int nb = (e.getsourceid() == id) ? e.gettargetid() : e.getsourceid();
            double w = e.getweight();
            if (perturb_ratio > 0.0) w *= noise(rng);
            csr.adjncy.push_back((idx_t)nb);
            csr.adjwgt.push_back((idx_t)max(1LL, llround(w * 1000.0)));
        }
        csr.xadj.push_back((idx_t)csr.adjncy.size());
    }
    csr.vwgt.assign(csr.nvtxs, 1);
    return csr;
}

static double compute_modularity(Graph& g, const vector<idx_t>& part, int num_parts)
{
    double intra_w = 0.0, total_w = 0.0;
    for (auto& [eid, e] : g.edgeset) {
        int s = e.getsourceid(), t = e.gettargetid();
        double w = e.getweight();
        total_w += w;
        if (s < (int)part.size() && t < (int)part.size() && part[s] == part[t])
            intra_w += w;
    }
    if (total_w < 1e-12) return 0.0;
    double Q = intra_w / total_w;

    vector<int> sz(num_parts, 0);
    for (auto v : part) if (v >= 0 && v < num_parts) sz[v]++;

    int empty_count = 0;
    for (int c : sz) if (c == 0) empty_count++;
    if (empty_count > 0) return Q - 1.0 * empty_count;

    double mean = (double)part.size() / num_parts;
    double var = 0.0;
    for (int c : sz) var += (c - mean) * (c - mean);
    double cv = sqrt(var / num_parts) / (mean + 1e-9);
    int max_sz = *max_element(sz.begin(), sz.end());
    double imbalance = (double)max_sz / (mean + 1e-9);

    double penalty;
    if (imbalance <= 1.5)      penalty = 0.05 * cv;
    else if (imbalance <= 2.0) penalty = 0.15 * cv;
    else                       penalty = 0.40 * cv;

    return Q - penalty;
}

static void rebalance_partition(vector<idx_t>& part, int num_parts, int nvtxs)
{
    vector<vector<int>> members(num_parts);
    for (int i = 0; i < nvtxs; i++)
        if (part[i] >= 0 && part[i] < num_parts) members[part[i]].push_back(i);
    for (int p = 0; p < num_parts; p++) {
        if (!members[p].empty()) continue;
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

static void force_balance_partition(Graph& g, vector<idx_t>& part, int num_parts, int nvtxs)
{
    double mean = (double)nvtxs / num_parts;
    // ceil(mean) + 허용 추가 1개: 예) mean=2.75이면 max=3+1=4
    int max_allowed = (int)ceil(mean) + 1;

    vector<vector<int>> members(num_parts);
    for (int i = 0; i < nvtxs; i++)
        if (part[i] >= 0 && part[i] < num_parts) members[part[i]].push_back(i);

    bool changed = true;
    int iters = 0;
    while (changed && iters < 20) {
        changed = false; iters++;
        for (int p = 0; p < num_parts; p++) {
            while ((int)members[p].size() > max_allowed) {
                int smallest = 0;
                for (int q = 1; q < num_parts; q++)
                    if (members[q].size() < members[smallest].size()) smallest = q;
                if (smallest == p) break;
                int best_node = -1;
                double best_gain = -1e18;
                for (int node : members[p]) {
                    double to_s = 0, to_p = 0;
                    if (g.nodeset.count(node))
                        for (int eid : g.nodeset.at(node).edges) {
                            if (!g.edgeset.count(eid)) continue;
                            const Edge& e = g.edgeset.at(eid);
                            int nb = (e.getsourceid() == node) ? e.gettargetid() : e.getsourceid();
                            if (nb < 0 || nb >= nvtxs) continue;
                            if (part[nb] == (idx_t)smallest) to_s += e.getweight();
                            else if (part[nb] == (idx_t)p)   to_p += e.getweight();
                        }
                    if (to_s - to_p > best_gain) { best_gain = to_s - to_p; best_node = node; }
                }
                if (best_node == -1) best_node = members[p].back();
                members[p].erase(find(members[p].begin(), members[p].end(), best_node));
                members[smallest].push_back(best_node);
                part[best_node] = (idx_t)smallest;
                changed = true;
            }
        }
    }
}

static bool run_metis_once(Graph& g, int num_parts, int seed,
                           double perturb_ratio, mt19937& rng,
                           int ufactor,
                           vector<idx_t>& part_out, idx_t& objval_out)
{
    MetisCSR csr = build_csr(g, perturb_ratio, rng);
    idx_t nvtxs = csr.nvtxs, ncon = 1, nparts = (idx_t)num_parts;
    part_out.assign(nvtxs, 0); objval_out = 0;

    idx_t opts[METIS_NOPTIONS];
    METIS_SetDefaultOptions(opts);
    opts[METIS_OPTION_SEED]    = (idx_t)((long long)seed * 1000003LL % 999983 + 1);
    opts[METIS_OPTION_NITER]   = 20;
    opts[METIS_OPTION_UFACTOR] = (idx_t)ufactor;

    int status = METIS_PartGraphKway(
        &nvtxs, &ncon,
        csr.xadj.data(), csr.adjncy.data(),
        csr.vwgt.data(), nullptr, csr.adjwgt.data(),
        &nparts, nullptr, nullptr,
        opts, &objval_out, part_out.data());

    if (status != METIS_OK) {
        MetisCSR csr2 = build_csr(g, 0.0, rng);
        idx_t nv2 = csr2.nvtxs;
        part_out.assign(nv2, 0); objval_out = 0;
        status = METIS_PartGraphKway(
            &nv2, &ncon,
            csr2.xadj.data(), csr2.adjncy.data(),
            csr2.vwgt.data(), nullptr, csr2.adjwgt.data(),
            &nparts, nullptr, nullptr,
            nullptr, &objval_out, part_out.data());
    }
    return (status == METIS_OK);
}

void Qcircuit::QMapper::run_metis_partition(Graph& interactiongraph, int num_parts)
{
    idx_t nvtxs = (idx_t)interactiongraph.nodeset.size();
    if (nvtxs <= (idx_t)num_parts) {
        partition_result.resize(nvtxs);
        for (idx_t i = 0; i < nvtxs; i++) partition_result[i] = i % num_parts;
        cout << "[ADAPTIVE-METIS] Trivial partition\n";
        return;
    }

    int num_edges = (int)interactiongraph.edgeset.size();
    double mean_sz = (double)nvtxs / num_parts;

    // ufactor: 크기별 — 소규모는 관대, 대규모는 엄격
    int ufactor;
    if      (nvtxs <= 12) ufactor = 1000;
    else if (nvtxs <= 20) ufactor = 500;
    else {
        int ceil_sz = (int)ceil(mean_sz);
        int max_sz_allowed = max(ceil_sz + 2, (int)ceil(mean_sz * 1.4));
        ufactor = max(150, (int)round(((double)max_sz_allowed / mean_sz - 1.0) * 1000));
        ufactor = min(ufactor, 300);
    }

    int MAX_TRIALS, PATIENCE, PURE_TRIALS;
    double MAX_PERTURB;
    if (num_edges <= 30) {
        MAX_TRIALS = 1; PATIENCE = 1; PURE_TRIALS = 1; MAX_PERTURB = 0.0;
    } else if (num_edges <= 100) {
        MAX_TRIALS = 50; PATIENCE = 12; PURE_TRIALS = 2; MAX_PERTURB = 0.40;
    } else {
        MAX_TRIALS = 100; PATIENCE = 18; PURE_TRIALS = 3; MAX_PERTURB = 0.45;
    }
    double IMPROVE_THRESH = 3e-4;

    double best_Q = -1e18;
    idx_t best_obj = (idx_t)INT_MAX;
    vector<idx_t> best_part;
    int no_improve = 0, trial = 0, success = 0;
    mt19937 rng(42);

    cout << "\n[ADAPTIVE-METIS] start: nvtxs=" << nvtxs
         << " edges=" << num_edges
         << " ufactor=" << ufactor
         << " max_trials=" << MAX_TRIALS << "\n";

    while (trial < MAX_TRIALS) {
        double perturb = 0.0;
        if (trial >= PURE_TRIALS) {
            double prog = (double)(trial - PURE_TRIALS) / max(1, MAX_TRIALS - PURE_TRIALS);
            perturb = MAX_PERTURB * prog;
        }
        int seed = trial * 137 + 42;
        vector<idx_t> part; idx_t objval;
        bool ok = run_metis_once(interactiongraph, num_parts, seed,
                                 perturb, rng, ufactor, part, objval);
        if (!ok) { trial++; continue; }
        success++;

        rebalance_partition(part, num_parts, (int)nvtxs);

        {
            vector<int> sz2(num_parts, 0);
            for (auto v : part) if (v >= 0 && v < num_parts) sz2[v]++;
            int cur_max = *max_element(sz2.begin(), sz2.end());
            if ((double)cur_max / (mean_sz + 1e-9) > 1.3)
                force_balance_partition(interactiongraph, part, num_parts, (int)nvtxs);
        }

        double Q = compute_modularity(interactiongraph, part, num_parts);

        if (best_part.empty() || Q > best_Q + IMPROVE_THRESH) {
            double prev_Q = best_Q;
            best_Q = Q; best_obj = objval; best_part = part; no_improve = 0;
            vector<int> sz(num_parts, 0);
            for (auto v : part) if (v >= 0 && v < num_parts) sz[v]++;
            cout << "  [trial " << setw(3) << trial << "] Q=" << fixed << setprecision(5) << Q
                 << " cut=" << objval << " perturb=" << setprecision(3) << perturb
                 << " sizes=[";
            for (int i = 0; i < num_parts; i++) cout << sz[i] << (i+1<num_parts?",":"");
            cout << "]";
            if (prev_Q > -1e17) cout << " Δ=+" << setprecision(5) << (Q - prev_Q);
            cout << "  *** BEST ***\n";
        } else {
            no_improve++;
            if (trial % 10 == 9)
                cout << "  [trial " << setw(3) << trial << "] Q=" << fixed << setprecision(5) << Q
                     << " best=" << best_Q << " no_imp=" << no_improve << "/" << PATIENCE << "\n";
        }
        if (no_improve >= PATIENCE) {
            cout << "  [ADAPTIVE-METIS] saturated at trial=" << trial << "\n";
            break;
        }
        trial++;
    }

    if (best_part.empty()) {
        cerr << "[ADAPTIVE-METIS] All trials failed, fallback\n";
        best_part.resize(nvtxs);
        for (idx_t i = 0; i < nvtxs; i++) best_part[i] = i % num_parts;
    }
    partition_result = std::move(best_part);

    vector<int> sz(num_parts, 0);
    for (auto v : partition_result) if (v >= 0 && v < num_parts) sz[v]++;
    int fmax = *max_element(sz.begin(), sz.end());
    cout << "[ADAPTIVE-METIS] DONE: best_Q=" << fixed << setprecision(5) << best_Q
         << " trials=" << trial+1 << "/" << MAX_TRIALS << " success=" << success << "\n"
         << "  partition sizes: [";
    for (int i = 0; i < num_parts; i++) cout << sz[i] << (i+1<num_parts?",":"");
    cout << "] imbalance=" << setprecision(2) << (double)fmax/(mean_sz+1e-9) << "\n";
}

void Qcircuit::QMapper::analyze_cross_partition_edges(
    Graph& g, const vector<idx_t>& part, int num_qpu)
{
    vector<vector<double>> cross(num_qpu, vector<double>(num_qpu, 0.0));
    for (const auto& [eid, e] : g.edgeset) {
        int s = e.getsourceid(), t = e.gettargetid();
        if (s < 0 || t < 0 || s >= (int)part.size() || t >= (int)part.size()) continue;
        int ps = part[s], pt = part[t];
        double w = e.getweight();
        if (ps == pt) continue;
        cross[ps][pt] += w; cross[pt][ps] += w;
    }
    cross_table = cross;
    double total_cross = 0.0;
    for (int i = 0; i < num_qpu; i++)
        for (int j = i+1; j < num_qpu; j++) total_cross += cross[i][j];
    cout << "\n========== Cross-Subgraph Weighted Sums ==========\n"
         << "  Total cross-QPU weight: " << fixed << setprecision(3) << total_cross << "\n";
    vector<tuple<int,int,double>> cl;
    for (int i = 0; i < num_qpu; i++)
        for (int j = i+1; j < num_qpu; j++) cl.push_back({i, j, cross[i][j]});
    sort(cl.begin(), cl.end(), [](auto& a, auto& b){ return get<2>(a) > get<2>(b); });
    for (auto& [i, j, val] : cl)
        if (val > 1e-9) cout << "  QPU" << i << "<->QPU" << j << ": " << setprecision(3) << val << "\n";
}

void Qcircuit::QMapper::reconstruct_info(Graph& interactiongraph, int num_parts)
{
    Subsets.clear(); Subsets.resize(num_parts);
    for (auto& [id, node] : interactiongraph.nodeset) {
        int p = partition_result[id];
        if (p >= 0 && p < num_parts) Subsets[p].nodeset[id] = node;
    }
    for (auto& [eid, e] : interactiongraph.edgeset) {
        int s = e.getsourceid(), t = e.gettargetid();
        int ps = partition_result[s], pt = partition_result[t];
        if (ps == pt && ps >= 0 && ps < num_parts) {
            Subsets[ps].edgeset[eid] = e;
            Subsets[ps].nodeset[s].edges.insert(eid);
            Subsets[ps].nodeset[t].edges.insert(eid);
        }
    }
    cout << "\n===== InteractionGraph Subsets Complete =====\n";
    inter_Sub_edge.clear();
    cout << "\n===== Building inter_Sub_edge (Cross-Sub Edges) =====\n";
    for (auto& [eid, e] : interactiongraph.edgeset) {
        int s = e.getsourceid(), t = e.gettargetid();
        int ps = partition_result[s], pt = partition_result[t];
        if (ps != pt) inter_Sub_edge.push_back({ ps, { pt, e } });
    }
}

void Qcircuit::QMapper::matching()
{
    cout << "\n===== MATCHING START =====\n";
    int num_sub = (int)Subsets.size();
    if (num_sub == 0) { cout << "No Subsets found.\n"; return; }
    matching_info.clear(); matching_info.resize(num_sub, { Graph(), -1 });
    for (int i = 0; i < num_sub; i++) {
        matching_info[i] = { Subsets[i], i };
        cout << "[MATCH] Sub" << i << " -> QPU" << i << "\n";
    }
    cout << "\n===== MATCHING DONE =====\n";
}