#include "circuit.h"
#include <vector>
#include<iostream>
#include <cmath>
using namespace std;
using namespace Qcircuit;
// 실수 가중치를 METIS용 정수 가중치로 변환해주는 함수

void Qcircuit::QMapper::analyze_cross_partition_edges(
    Graph& g,
    const vector<idx_t>& part,
    int num_qpu)
{
    // QPU 간 가중치 합 누적용 행렬

   
    vector<vector<double>> cross(num_qpu, vector<double>(num_qpu, 0.0));
    // interactionGraph의 모든 edge 순회
    for (const auto& [eid, e] : g.edgeset)
    {
        int s = e.getsourceid();
        int t = e.gettargetid();

        // part 인덱스 범위 체크
        if (s < 0 || t < 0 ||
            s >= static_cast<int>(part.size()) ||
            t >= static_cast<int>(part.size()))
            continue;

        int ps = part[s];
        int pt = part[t];

        double w = e.getweight();

        if (ps == pt) continue; // 같은 파티션이면 skip

        cross[ps][pt] += w;
        cross[pt][ps] += w; // 무방향 그래프라 대칭 누적
    }
    cross_table=cross;

    // ===== 출력 =====
    cout << "\n========== Cross-QPU Weighted Edge Matrix ==========\n";
    cout << "      ";
    for (int j = 0; j < num_qpu; j++)
        cout << "  Sub" << j << " ";
    cout << "\n";

    for (int i = 0; i < num_qpu; i++)
    {
        cout << "Sub" << i << " ";
        for (int j = 0; j < num_qpu; j++)
            cout << setw(8) << cross[i][j] << " ";
        cout << "\n";
    }

    // 정렬용 리스트
    vector<tuple<int,int,double>> cross_list;
    for (int i = 0; i < num_qpu; i++)
        for (int j = i + 1; j < num_qpu; j++)
            cross_list.push_back({i, j, cross[i][j]});

    sort(cross_list.begin(), cross_list.end(),
         [](auto &a, auto &b){ return get<2>(a) > get<2>(b); });

    cout << "\n========== Cross-Subgraph Weighted Sums (Sorted) ==========\n";
    for (auto &tup : cross_list)
    {
        int i, j;
        double val;
        tie(i, j, val) = tup;
        cout << "Sub" << i << " - Sub" << j
             << " : " << val << " Total weight\n";
    }
}


void Qcircuit::QMapper::run_metis_partition(Graph& interactiongraph, int num_parts)
{
    // -----------------------------
    // 기본 메타 정보
    // -----------------------------
    idx_t nvtxs_1 = interactiongraph.nodeset.size();  // 노드 수
    idx_t ncon_1  = 1;                                // constraint 수

    // -----------------------------
    // CSR 포맷용 배열 준비
    // -----------------------------
    vector<idx_t> xadj_1;
    vector<idx_t> adjncy_1;
    vector<idx_t> adjwgt_1;

    xadj_1.reserve(nvtxs_1 + 1);
    xadj_1.push_back(0);

    for (auto& [id, node] : interactiongraph.nodeset)
    {
        for (int edge_id : node.edges)
        {
            Edge& e = interactiongraph.edgeset[edge_id];
            int neighbor = (e.getsourceid() == id) ? e.gettargetid()
                                                   : e.getsourceid();

            adjncy_1.push_back(neighbor);
            adjwgt_1.push_back(e.getweight());
        }
        xadj_1.push_back(adjncy_1.size());
    }

    // -----------------------------
    // 노드 가중치 (없으면 기본 1)
    // -----------------------------
    vector<idx_t> vwgt_1(nvtxs_1, 1);

    // -----------------------------
    // METIS 결과 저장 변수들
    // -----------------------------
    idx_t nparts_1 = num_parts;
    idx_t objval_1 = 0;
    vector<idx_t> part(nvtxs_1, 0);

    // -----------------------------
    // METIS 실행
    // -----------------------------
    int status = METIS_PartGraphKway(
        &nvtxs_1,
        &ncon_1,
        xadj_1.data(),
        adjncy_1.data(),
        vwgt_1.data(),
        /* vsize */ NULL,
        adjwgt_1.data(),
        &nparts_1,
        /* tpwgts */ NULL,
        /* ubvec  */ NULL,
        /* options*/ NULL,
        &objval_1,
        part.data()
    );

    // -----------------------------
    // 결과 처리
    // -----------------------------
    if (status == METIS_OK)
    {
        cout << "\n[METIS] Partitioning complete\n";
        cout << "Total Sum of Weight: " << objval_1 << "\n";

        for (int i = 0; i < nvtxs_1; i++)
            cout << "Node " << setw(2) << i << " → Part " << part[i] << "\n";

        partition_result = std::move(part);
    }
    else
    {
        cerr << "[METIS] Partitioning failed.\n";
    }
 
}



void Qcircuit::QMapper::reconstruct_info(Graph& interactiongraph, int num_parts)
{   Subsets.clear();
    Subsets.resize(num_parts);
    // 1) 노드 분배
for (auto &[id, node] : interactiongraph.nodeset)
{
    int p = partition_result[id]; // 그냥 id 를 part 인덱스로 사용
    if (p >= 0 && p < num_parts)
    {
        Subsets[p].nodeset[id] = node;   // 노드 그대로 복사
    }
}

// 2) 엣지 분배 (양 노드가 같은 파티션에 있을 때만)
for (auto &[eid, e] : interactiongraph.edgeset)
{
    int s = e.getsourceid();
    int t = e.gettargetid();

    int ps = partition_result[s];
    int pt = partition_result[t];

    if (ps == pt && ps >= 0 && ps < num_parts)
    {
        // edge 넣기
        Subsets[ps].edgeset[eid] = e;

        // 해당 subgraph에서 node-edge 관계 구축
        Subsets[ps].nodeset[s].edges.insert(eid);
        Subsets[ps].nodeset[t].edges.insert(eid);
    }
}

cout << "\n===== InteractionGraph Subsets Complete =====\n";
for (int p = 0; p < num_parts; p++)
{
    cout << "[Subset " << p << "]\n";
    cout << "Nodes: " << Subsets[p].nodeset.size() 
         << "  Edges: " << Subsets[p].edgeset.size() << "\n\n";
}
inter_Sub_edge.clear();

cout << "\n===== Building inter_Sub_edge (Cross-Sub Edges) =====\n";

// interactionGraph의 모든 edge 순회
for (auto &[eid, e] : interactiongraph.edgeset)
{
    int s  = e.getsourceid();
    int t  = e.gettargetid();

    int ps = partition_result[s];   
    int pt = partition_result[t];   

    // 서로 다른 QPU면 inter-QPU edge
    if (ps != pt)
    {
        
        inter_Sub_edge.push_back({ ps, { pt, e } });

        cout << "[INTER-Sub] EdgeID " << eid 
             << " : q" << s << "(Sub" << ps << ")"
             << " -- q" << t << "(Sub" << pt << ")"
             << "  weight=" << e.getweight()
             << endl;
    }
}

cout << "Total inter-Subgraph edges: " << inter_Sub_edge.size() << "\n\n";

// ===== Optional: 구조 점검 =====
cout << "===== inter_Subgraph Content =====\n";
for (auto &info : inter_Sub_edge)
{
    int srcSub = info.first;
    int dstSub = info.second.first;
    const Edge &edgeObj = info.second.second;

    cout << " Sub" << srcSub 
         << " → Sub" << dstSub 
         << " | edge (" << edgeObj.getsourceid()
         << " -- " << edgeObj.gettargetid() 
         << "), weight=" << edgeObj.getweight()
         << ", EdgeID=" << edgeObj.getid()
         << "\n";
}
cout << "===================================\n";

}

void Qcircuit::QMapper::matching()
{
     cout << "\n===== MATCHING START =====\n";

    int num_sub = Subsets.size();
    if (num_sub == 0)
    {
        cout << "No Subsets found.\n";
        return;
    }

    // ----------------------------
    // cross_table 기반 Sub 중요도 계산
    // ----------------------------
    vector<pair<int,double>> importance; 
    importance.resize(num_sub);

    for (int i = 0; i < num_sub; i++)
    {
        importance[i].first = i;
        importance[i].second = 0.0;

        for (int j = 0; j < num_sub; j++)
        {
            if (i == j) continue;
            importance[i].second += cross_table[i][j];  // 중요도 계산
        }
    }

    // 중요도 높은 순으로 정렬
    sort(importance.begin(), importance.end(),
        [](auto &a, auto &b){
            return a.second > b.second;
        });

    // ----------------------------
    // degree 높은 QPU 순서대로 배치
    // ----------------------------
    matching_info.clear();
    matching_info.resize(num_sub, {Graph(), -1});

    vector<bool> used(num_sub, false);

    for (auto &[subIdx, imp] : importance)
    {
        if (qpu_degree_graph.degree_queue.empty()) break;

        auto [qpu, deg] = qpu_degree_graph.degree_queue.front();
        qpu_degree_graph.degree_queue.pop();

        matching_info[subIdx] = { Subsets[subIdx], qpu };
        used[subIdx] = true;

        cout << "[MATCH] Sub" << subIdx 
             << " (sumWeight=" << imp 
             << ") -> QPU" << qpu << "\n";
    }

    // ----------------------------
    // 결과 출력
    // ----------------------------
    cout << "\n===== MATCHING DONE =====\n";
    for (int i = 0; i < num_sub; i++)
    {
        cout << "Sub" << i << " → QPU" << matching_info[i].second << "\n";
    }
}