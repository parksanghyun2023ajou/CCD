#include "circuit.h"

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::generate_multi_qpu(int num_qpu, int num_qubits, bool fully_connected_qubit, bool fully_connected_qpu, bool buffer_insertion)
{
  cout << "Generating Multi-QPU Architecture: " << num_qpu
<< " Compute QPUs + 1 Buffer QPU, each with " << num_qubits << " qubits\n";

multi_qpu_graph = Graph();

int a = static_cast<int>(sqrt(num_qubits));
if (a * a != num_qubits) {
    cerr << "Error: num_qubits must be a perfect square.\n";
    exit(1);
}

int total_qpu_count = num_qpu + 1;
int buffer_qpu_idx = num_qpu; // 버퍼 코어는 맨 마지막 인덱스 사용

int total_nodes = total_qpu_count * num_qubits;

for (int i = 0; i < total_nodes; i++) {
    multi_qpu_graph.addnode();
}

// 1. 코어 내부 토폴로지 생성

for (int qpu_idx = 0; qpu_idx < total_qpu_count; qpu_idx++) {
    int base = qpu_idx * num_qubits;

    // 버퍼 코어일 경우
    if (qpu_idx == buffer_qpu_idx) {
        // Fully connected across different rows
        for (int r1 = 0; r1 < a; r1++) {
            for (int c1 = 0; c1 < a; c1++) {
                int u = base + r1 * a + c1;

                // 중복 연결을 피하기 위해 r2는 r1보다 큰 행만 탐색
                for (int r2 = r1 + 1; r2 < a; r2++) {
                    for (int c2 = 0; c2 < a; c2++) {
                        int v = base + r2 * a + c2; // 다른 행에 있는 큐비트
                        multi_qpu_graph.addedge(u, v, 1.0);
                    }
                }
            }
        }
    }
    // 일반 연산 코어일 경우
    else {
        for (int r = 0; r < a; r++) {
            for (int c = 0; c < a; c++) {
                int id = base + r * a + c;

                // 우측 연결
                if (c + 1 < a)
                    multi_qpu_graph.addedge(id, id + 1, 1.0);

                // 아래 연결
                if (r + 1 < a)
                    multi_qpu_graph.addedge(id, id + a, 1.0);

                // 대각선 연결 (옵션)
                if (fully_connected_qubit && r + 1 < a && c + 1 < a) {
                    multi_qpu_graph.addedge(id, id + a + 1, 1.0);
                    multi_qpu_graph.addedge(id + 1, id + a, 1.0);
                }
            }
        }
    }
}

// 2. Load / Unload 통신망 구축 (Inter-QPU Edges)

double inter_weight = buffer_insertion ? 3.0 : 5.0;
int buffer_base = buffer_qpu_idx * num_qubits;

for (int i = 0; i < num_qpu; i++) {

    // 연산 코어 개수가 버퍼의 행 개수(a)보다 많으면 연결 불가
    if (i >= a) {
        cerr << "\n[Warning] Not enough rows in Buffer QPU to connect QPU " << i << "!\n";
        break;
    }

    int qpu_base = i * num_qubits;

    // QPUi와 상단 row와 버퍼 코어의 할당된 행 연결
    for (int c1 = 0; c1 < a; c1++) {
        int qpu_node = qpu_base + c1; // 연산 코어의 최상단 노드들

        for (int c2 = 0; c2 < a; c2++) {
            int buf_node = buffer_base + (i * a) + c2; // 버퍼 코어에서 나에게 할당된 행의 노드들

            multi_qpu_graph.addedge(qpu_node, buf_node, inter_weight);
        }
    }
}


// 3. QPU Degree Graph

cout << "\nSkipping complex QPU Degree Graph (Star Topology applied)\n";
qpu_degree_graph = Graph();

// 연산 코어들만 순서대로 
qpu_degree_graph.degree_queue = queue<pair<int, int>>();
for (int i = 0; i < num_qpu; i++) {
    qpu_degree_graph.degree_queue.push({ i, 1 });
}

/*//아래는 검증용 검증완료
cout << "\n=== Multi-QPU Interconnect Edges ===\n";

for (const auto& [eid, e] : multi_qpu_graph.edgeset)
{
    int s = e.getsourceid();
    int t = e.gettargetid();

    int qpu_s = s / num_qubits;
    int qpu_t = t / num_qubits;

    int local_s = s % num_qubits;
    int local_t = t % num_qubits;

    cout << "Edge " << eid << " : ";

    if (qpu_s == buffer_qpu_idx)
        cout << "[Buffer QPU " << qpu_s << "]";
    else
        cout << "[QPU " << qpu_s << "]";

    cout << ".q" << local_s << " <--> ";

    if (qpu_t == buffer_qpu_idx)
        cout << "[Buffer QPU " << qpu_t << "]";
    else
        cout << "[QPU " << qpu_t << "]";

    cout << ".q" << local_t
        << "  (weight=" << e.getweight() << ")\n";
}*/


}