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



// ========================================
// 5) initial_mapping() 본 함수 (SAFE)
// ========================================
void Qcircuit::QMapper::initial_mapping(int num_qpu, int num_qubit)
{
    cout << "initial_mapping\n";

    circuit_processing(num_qpu);
    this->positions = num_qubit;

    layout_L.clear();
    qubit_Q.clear();
    layout_L.resize(num_qpu);
    qubit_Q.resize(num_qpu);

    for (int sid = 0; sid < (int)matching_info.size(); sid++)
    {
        Graph& subG = matching_info[sid].first;
        int qpu = matching_info[sid].second;

        cout << "\n-------------------------------------------------------------\n";
        cout << "[Subgraph " << sid << " → QPU " << qpu << "]\n";

        if (subG.nodeset.empty()) {
            cout << "(empty subgraph)\n";
            continue;
        }

        cout << "Global logical node IDs: ";
        for (auto &p : subG.nodeset)
            cout << p.first << " ";
        cout << "\n";

        // -------------------------
        // logical ordering
        // -------------------------
        int logical_center = pick_logical_center(subG);
        vector<int> logical_order;
        bfs_order_logical(subG, logical_center, logical_order);

        cout << "Logical center: " << logical_center << "\n";

        // -------------------------
        // physical ordering
        // -------------------------
        int phys_center = pick_physical_center(multi_qpu_graph, qpu, positions);
        vector<int> physical_order;
        bfs_order_physical(multi_qpu_graph, qpu, positions, phys_center, physical_order);

        cout << "Physical center (QPU " << qpu << "): " << phys_center << "\n";

        // -------------------------
        // mapping
        // -------------------------
        int max_map = min((int)logical_order.size(),
                          (int)physical_order.size());

        cout << "[MAP]\n";
        for (int i = 0; i < max_map; i++)
        {
            int logical = logical_order[i];
            int physical = physical_order[i];

            layout_L[qpu][logical] = physical;
            qubit_Q[qpu][physical] = logical;

            cout << " logical " << logical
                 << "  →  physical " << physical << "\n";
        }
    }

     cout << "\n===================== initial_mapping END =====================\n";
}

