// ddd
// new_generate_qpu.cpp

#include "circuit.h"
#include <cmath> 

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::generate_multi_qpu(int num_qpu, int num_qubits, bool fully_connected_qubit, bool buffer_insertion)
{
    cout << "Generating Multi-QPU Architecture: " << num_qpu
         << " QPUs, each with " << num_qubits << " qubits\n";

    multi_qpu_graph = Graph();

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    if (n * (n) != num_qubits) {
        cerr << "Error: num_qubits must be of the form (n+1) * n\n";
        exit(1);
    }

    int rows = n + 1;
    int cols = n;

    int total_nodes = num_qpu * (num_qubits+n); 

    for (int i = 0; i < total_nodes; i++) {
        multi_qpu_graph.addnode();
    }

    // 1. 코어 내부 토폴로지
    for (int qpu_idx = 0; qpu_idx < num_qpu; qpu_idx++) {
        int base = qpu_idx * (num_qubits+n); 

        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                int id = base + r * cols + c;

                // 우측 연결
                if (c + 1 < cols&&r!=0)
                    multi_qpu_graph.addedge(id, id+1, 1.0);

                // 아래 연결 (연산 노드와 버퍼 노드 간의 길)
                if (r + 1 < rows)
                    multi_qpu_graph.addedge(id, id + cols, 1.0);

                // 대각선 연결 (옵션)
                if (fully_connected_qubit && r + 1 < rows && c + 1 < cols) {
                    multi_qpu_graph.addedge(id, id + cols + 1, 1.0);
                    multi_qpu_graph.addedge(id + 1, id  + cols, 1.0);
                }
            }
        }
    }

    // 2. QPU Degree Graph
    cout << "\nSkipping complex QPU Degree Graph (Using simple Star Topology meta-info)\n";
    qpu_degree_graph = Graph(); 

    qpu_degree_graph.degree_queue = queue<pair<int, int>>();
    for (int i = 0; i < num_qpu; i++) {
        qpu_degree_graph.degree_queue.push({i, 1}); 
    }
    cout << "\n=== Multi-QPU Interconnect Edges ===\n";

for (const auto& [eid, e] : multi_qpu_graph.edgeset)
{
    int s = e.getsourceid();
    int t = e.gettargetid();

    int qpu_s = s / (num_qubits+n);
    int qpu_t = t / (num_qubits+n);

    int local_s = s % (num_qubits+n);
    int local_t = t % (num_qubits+n);

    cout << "Edge " << eid << " : ";
        cout << "[QPU " << qpu_s << "]";

    cout << ".q" << local_s << " <--> ";
        cout << "[QPU " << qpu_t << "]";

    cout << ".q" << local_t
        << "  (weight=" << e.getweight() << ")\n";
}
}