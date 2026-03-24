#include "circuit.h"

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::generate_multi_qpu(int num_qpu, int num_qubits, bool fully_connected_qubit, bool fully_connected_qpu, bool buffer_insertion)
{
    cout << "Generating Multi-QPU Architecture: " << num_qpu
         << " QPUs, each with " << num_qubits << " qubits\n";

    multi_qpu_graph = Graph(); 

    int a = static_cast<int>(sqrt(num_qubits));
    if (a * a != num_qubits) {
        cerr << "Error: num_qubits must be a perfect square.\n";
        exit(1);
    }

    int total_nodes = num_qpu * num_qubits;

    for (int i = 0; i < total_nodes; i++) {
        multi_qpu_graph.addnode();
    }

    for (int qpu_idx = 0; qpu_idx < num_qpu; qpu_idx++) {
        int base = qpu_idx * num_qubits;

        for (int r = 0; r < a; r++) {
            for (int c = 0; c < a; c++) {
                int id = base + r * a + c;

                if (c + 1 < a)
                    multi_qpu_graph.addedge(id, id + 1, 1.0);

                if (r + 1 < a)
                    multi_qpu_graph.addedge(id, id + a, 1.0);

                if (fully_connected_qubit && r + 1 < a && c + 1 < a) {
                    multi_qpu_graph.addedge(id, id + a + 1, 1.0);
                    multi_qpu_graph.addedge(id + 1, id + a, 1.0);
                }
            }
        }
    }

    double inter_weight = buffer_insertion ? 3.0 : 5.0;
    int interconnect_qubits = a;


for (int i = 0; i < num_qpu; i++) {

    int basis;

    if (fully_connected_qpu) {
        // Ring
        basis = (i + 1) % num_qpu;

        // 2개 QPU 예외 → 한 번만 연결
        if (num_qpu == 2 && basis == 0)
            continue;

    } else {
        // Linear chain
        if (i == num_qpu - 1)
            break;

        basis = i + 1;
    }

    int base_i = i * num_qubits;
    int base_j = basis * num_qubits;

    for (int r = 0; r < interconnect_qubits; r++) {

        int from = base_i + r;
        int to   = base_j + r;

        multi_qpu_graph.addedge(from, to, inter_weight);

        // 검증 완료, 아래는 검증용 코드임
        /*int qpu1   = from / num_qubits;
        int qubit1 = from % num_qubits;

        int qpu2   = to / num_qubits;
        int qubit2 = to % num_qubits;

        cout
            << "<" << qpu1 << "," << qubit1 << ">"
            << " <=> "
            << "<" << qpu2 << "," << qubit2 << ">"
            << endl;*/
    }

    }

    cout << "\nBuilding QPU Degree Graph\n";

    qpu_degree_graph = Graph(); 

    for (int i = 0; i < num_qpu; ++i) {
        qpu_degree_graph.addnode(); 
    }

    set<pair<int, int>> qpus_connections;

    for (const auto& [edge_id, edge] : multi_qpu_graph.edgeset) {
        int source_qubit_id = edge.getsourceid();
        int target_qubit_id = edge.gettargetid();

        int qpu_of_source = source_qubit_id / num_qubits;
        int qpu_of_target = target_qubit_id / num_qubits;

        if (qpu_of_source != qpu_of_target) { 
            int qpu_i = min(qpu_of_source, qpu_of_target);
            int qpu_j = max(qpu_of_source, qpu_of_target);
            qpus_connections.insert({qpu_i, qpu_j});
        }
    }

    for (const auto& conn : qpus_connections) {
        qpu_degree_graph.addedge(conn.first, conn.second, 1.0);
    }

    qpu_degree_graph.degree_info.clear();
    qpu_degree_graph.degree_info.reserve(num_qpu);

    for (int i = 0; i < num_qpu; i++)
    {
        int degree = qpu_degree_graph.nodeset[i].edges.size();
        qpu_degree_graph.degree_info.emplace_back(i, degree);
    }

    vector<pair<int,int>> degree_info_unsorted = qpu_degree_graph.degree_info;

    sort(qpu_degree_graph.degree_info.begin(), qpu_degree_graph.degree_info.end(),
         [](const pair<int,int>& a, const pair<int,int>& b){
             return a.second > b.second;
         });

    qpu_degree_graph.degree_queue = queue<pair<int, int>>(); 
    
    for (const auto& p : qpu_degree_graph.degree_info)
    {
        qpu_degree_graph.degree_queue.push(p); 
    }

    cout << "\n=== Original degree_info (stored in Graph) ===\n";
    for (const auto& p : degree_info_unsorted)
    {
        cout << "QPU " << p.first
             << "  Degree = " << p.second << "\n";
    }

    cout << "\n=== degree_queue (FIFO, stored in Graph) ===\n";
    {
        queue<pair<int,int>> tmp = qpu_degree_graph.degree_queue; 

        while (!tmp.empty())
        {
            auto [id, deg] = tmp.front();
            cout << "QPU " << id
                 << "  Degree = " << deg << "\n";
            tmp.pop();
        }
    }
    multi_qpu_graph.build_dist_table();
}