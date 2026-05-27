// initial_mapping.cpp

#include "circuit.h"
#include <queue>
#include <unordered_map>
#include <algorithm>
#include <climits>
#include <iostream>
// dododo

using namespace std;
using namespace Qcircuit;

namespace {

    // ======================
    // 1) Logical center pick
    // ======================
    int pick_logical_center(const Graph& g)
    {
        int best_id = -1;
        int best_deg = -1;

        for (auto& kv : g.nodeset) {
            int id = kv.first;
            int deg = kv.second.edges.size();
            if (deg > best_deg) {
                best_deg = deg;
                best_id = id;
            }
        }
        if (best_id == -1 && !g.nodeset.empty())
            best_id = g.nodeset.begin()->first;
        return best_id;
    }

    // ========================================
    // 2) BFS order for logical (SAFE VERSION)
    // ========================================
    void bfs_order_logical(const Graph& g,
                           int center,
                           vector<int>& order_out)
    {
        order_out.clear();
        if (g.nodeset.empty())
            return;

        unordered_map<int,int> dist;
        unordered_map<int,int> degree;

        for (auto& kv : g.nodeset)
            degree[kv.first] = kv.second.edges.size();

        queue<int> q;
        dist[center] = 0;
        q.push(center);

        while (!q.empty())
        {
            int u = q.front(); q.pop();

            const Node& n = g.nodeset.at(u);

            for (int eid : n.edges)
            {
                // 🔥 Subgraph에 존재하는 edge만 사용
                if (!g.edgeset.count(eid))
                    continue;

                const Edge& e = g.edgeset.at(eid);
                int s = e.source->getid();
                int t = e.target->getid();
                int v = (s == u ? t : s);

                if (!g.nodeset.count(v)) 
                    continue;
                if (dist.count(v)) 
                    continue;

                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }

        vector<tuple<int,int,int>> tmp;
        for (auto& kv : g.nodeset) {
            int id = kv.first;
            int d = (dist.count(id) ? dist[id] : INT_MAX/2);
            tmp.emplace_back(d, -degree[id], id);
        }

        sort(tmp.begin(), tmp.end(),
            [](auto& a, auto& b){
                if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
                if (get<1>(a) != get<1>(b)) return get<1>(a) < get<1>(b);
                return get<2>(a) < get<2>(b);
            });

        for (auto& t : tmp)
            order_out.push_back(get<2>(t));
    }

    // ========================================
    // 3) pick physical center in one QPU
    // ========================================
    int pick_physical_center(const Graph& multi_qpu_graph,
                             int qpu,
                             int positions)
    {
        int base = qpu * positions;
        int end  = base + positions;

        int best_id = base;
        int best_deg = -1;

        for (int p = base; p < end; p++)
        {
            if (!multi_qpu_graph.nodeset.count(p))
                continue;

            const Node& node = multi_qpu_graph.nodeset.at(p);
            int deg = 0;

            for (int eid : node.edges)
            {
                if (!multi_qpu_graph.edgeset.count(eid))
                    continue;

                const Edge& e = multi_qpu_graph.edgeset.at(eid);
                int s = e.source->getid();
                int t = e.target->getid();
                int v = (s == p ? t : s);

                if (v >= base && v < end)
                    deg++;
            }

            if (deg > best_deg) {
                best_deg = deg;
                best_id = p;
            }
        }
        return best_id;
    }

    // ========================================
    // 4) BFS order for physical (SAFE VERSION)
    // ========================================
    void bfs_order_physical(const Graph& mg,
                             int qpu,
                             int positions,
                             int center_phys,
                             vector<int>& order_out)
    {
        order_out.clear();
        if (mg.nodeset.empty())
            return;

        int base = qpu * positions;
        int end  = base + positions;

        unordered_map<int,int> dist;
        unordered_map<int,int> degree;

        for (int p = base; p < end; p++) {
            if (mg.nodeset.count(p))
                degree[p] = mg.nodeset.at(p).edges.size();
        }

        queue<int> q;
        dist[center_phys] = 0;
        q.push(center_phys);

        while (!q.empty()) {
            int u = q.front(); q.pop();

            if (!mg.nodeset.count(u)) continue;
            const Node& n = mg.nodeset.at(u);

            for (int eid : n.edges)
            {
                if (!mg.edgeset.count(eid))
                    continue;

                const Edge& e = mg.edgeset.at(eid);
                int s = e.source->getid();
                int t = e.target->getid();
                int v = (s == u ? t : s);

                if (v < base || v >= end)
                    continue;
                if (dist.count(v))
                    continue;

                dist[v] = dist[u] + 1;
                q.push(v);
            }
        }

        vector<tuple<int,int,int>> tmp;
        for (int p = base; p < end; p++) {
            if (!mg.nodeset.count(p)) continue;
            int d = (dist.count(p) ? dist[p] : INT_MAX/2);
            tmp.emplace_back(d, -degree[p], p);
        }

        sort(tmp.begin(), tmp.end(),
            [](auto& a, auto& b){
                if (get<0>(a) != get<0>(b)) return get<0>(a) < get<0>(b);
                if (get<1>(a) != get<1>(b)) return get<1>(a) < get<1>(b);
                return get<2>(a) < get<2>(b);
            });

        for (auto& t : tmp)
            order_out.push_back(get<2>(t));
    }

} // namespace



void Qcircuit::QMapper::initial_mapping(int num_qpu, int num_qubit)
{
    cout << "initial_mapping\n";

    circuit_processing(num_qpu);
    this->positions = num_qubit;

    layout_L.clear(); layout_L.resize(num_qpu);
    qubit_Q.clear();  qubit_Q.resize(num_qpu);
    for (int q = 0; q < num_qpu; q++) {
        layout_L[q].clear();
        qubit_Q[q].clear();
    }

    int n = static_cast<int>(std::floor(std::sqrt(num_qubit)));

    // =========================================================
    // STEP 1: logical → QPU 매핑 구성 (METIS 결과 기반)
    // =========================================================
    unordered_map<int,int> logical_to_qpu;
    for (int sid = 0; sid < (int)matching_info.size(); sid++) {
        Graph& subG = matching_info[sid].first;
        int qpu     = matching_info[sid].second;
        for (auto& kv : subG.nodeset)
            logical_to_qpu[kv.first] = qpu;
    }

    // =========================================================
    // STEP 2: inter_degree 계산
    //   Dgraph(전체 회로) 순회 → cross-QPU 게이트면 양쪽 누적
    // =========================================================
    vector<int> inter_degree(nqubits, 0);
    for (auto& gate : Dgraph.nodeset) {
        int c = gate.control;
        int t = gate.target;
        if (c < 0 || t < 0) continue;
        if (!logical_to_qpu.count(c) || !logical_to_qpu.count(t)) continue;
        if (logical_to_qpu[c] != logical_to_qpu[t]) {
            inter_degree[c]++;
            inter_degree[t]++;
        }
    }

    // =========================================================
    // STEP 3: 서브그래프별 배치
    // =========================================================
    for (int sid = 0; sid < (int)matching_info.size(); sid++)
    {
        Graph& subG = matching_info[sid].first;
        int    qpu  = matching_info[sid].second;
        if (subG.nodeset.empty()) continue;

        // 물리 슬롯: BFS 순서로 전체 정렬 후 top row / inner 분리
        int phys_center = pick_physical_center(multi_qpu_graph, qpu, num_qubit);
        vector<int> physical_order;
        bfs_order_physical(multi_qpu_graph, qpu, num_qubit,
                           phys_center, physical_order);

        vector<int> top_row_slots, inner_slots;
        for (int p : physical_order) {
            int row = (p % num_qubit) / n;
            if (row == 0) top_row_slots.push_back(p);
            else          inner_slots.push_back(p);
        }

        // 논리 큐비트 목록
        vector<int> logical_nodes;
        for (auto& kv : subG.nodeset)
            logical_nodes.push_back(kv.first);

        // inter_degree 내림차순 정렬
        sort(logical_nodes.begin(), logical_nodes.end(),
             [&](int a, int b){ return inter_degree[a] > inter_degree[b]; });

        // Top-3 선택 (inter_degree > 0, top row 슬롯 수 이내)
        int top_limit = min({3,
                             (int)top_row_slots.size(),
                             (int)logical_nodes.size()});

        vector<int> top3, rest;
        for (int i = 0; i < (int)logical_nodes.size(); i++) {
            if ((int)top3.size() < top_limit
                && inter_degree[logical_nodes[i]] > 0)
                top3.push_back(logical_nodes[i]);
            else
                rest.push_back(logical_nodes[i]);
        }

        // Top-3 → top row 배치
        for (int i = 0; i < (int)top3.size(); i++) {
            int logical  = top3[i];
            int physical = top_row_slots[i];
            layout_L[qpu][logical] = physical;
            qubit_Q[qpu][physical] = logical;
        }

        // =====================================================
        // STEP 4: BFS center 선택
        //   top3와 subG 에지 weight(게이트 빈도) 합산 최대 큐비트
        // =====================================================
        unordered_set<int> top3_set(top3.begin(), top3.end());

        int    bfs_center = -1;
        double best_score = -1.0;

        for (int q : rest) {
            if (!subG.nodeset.count(q)) continue;
            double score = 0.0;
            for (int eid : subG.nodeset.at(q).edges) {
                if (!subG.edgeset.count(eid)) continue;
                const Edge& e = subG.edgeset.at(eid);
                int other = (e.source->getid() == q)
                            ? e.target->getid()
                            : e.source->getid();
                if (top3_set.count(other))
                    score += e.getweight();   // 게이트 빈도(weight) 합산
            }
            // 동점이면 전체 inter_degree로 보조 비교
            if (score > best_score ||
                (score == best_score && bfs_center != -1 &&
                 inter_degree[q] > inter_degree[bfs_center]))
            {
                best_score = score;
                bfs_center = q;
            }
        }

        // fallback: top3와 연결 없으면 기존 방식
        if (bfs_center == -1)
            bfs_center = pick_logical_center(subG);

        // =====================================================
        // STEP 5: 나머지 → BFS 순서로 inner 슬롯 채우기
        // =====================================================
        if (!rest.empty() && bfs_center != -1) {
            vector<int> logical_order;
            bfs_order_logical(subG, bfs_center, logical_order);

            // inner 먼저, 남은 top row 나중
            vector<int> remaining_phys;
            for (int p : inner_slots)
                remaining_phys.push_back(p);
            for (int i = (int)top3.size(); i < (int)top_row_slots.size(); i++)
                remaining_phys.push_back(top_row_slots[i]);

            int pi = 0;
            for (int logical : logical_order) {
                if (top3_set.count(logical)) continue;  // 이미 배치됨
                if (pi >= (int)remaining_phys.size()) break;
                layout_L[qpu][logical] = remaining_phys[pi];
                qubit_Q[qpu][remaining_phys[pi]] = logical;
                pi++;
            }
        }
    }

    cout << "\n===================== initial_mapping END =====================\n";
}