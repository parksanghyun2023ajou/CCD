//For meet 수정본
//main_mapping.cpp
#include "circuit.h"

using namespace std;
using namespace Qcircuit;

#define Dlist_all_mode 0 


//////////////////////////// 추가한 함수 ////////////////////////
double Qcircuit::QMapper::mapping_machine(bool is_inter_gate, const pair<int, int> SWAP_pair, Circuit& dgraph, int gateid) {

    // 가상 스왑 진행할 두 물큐의 인덱스
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;
    int qpu_idx1 = Q1 / num_qubits; // Q1의 노드 id가 그대로 들어간다고 가정하고 qpu 인덱스 구현
    int qpu_idx2 = Q2 / num_qubits; //Q2의 노드 id가 그대로 들어간다고 가정하고 qpu 인덱스 구현

    // 해당 물큐에 올라가 있는 논큐 인덱스
    int q1 = qubit_Q[qpu_idx1][Q1];
    int q2 = qubit_Q[qpu_idx2][Q2];

    // 타깃 게이트의 타깃, 제어 큐비트 인덱스
    int control = dgraph.nodeset[gateid].control;
    int target = dgraph.nodeset[gateid].target;

    double final_cost = 0.0;

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    double WEIGHT_SWAP = 3.0;           // 코어 내부에서 1칸 SWAP 시 페널티
    double WEIGHT_VIRTUAL_BUFFER = 10.0; // 가상 버퍼에 올리는 기본 통신 페널티

    if (is_inter_gate)  //Inter의 경우
    {
        // Cost 계산: (최상단까지 가기위한 Swap 횟수 * SWAP 비용) + 가상 버퍼 업로드 비용

        // 그 논큐가 타코어와 연산해야하는 제어 또는 타깃인 경우
        if (q1 == control || q1 == target) {

            // 스왑 전: 현재 Q1 위치에서 최상단 Row까지 가기 위해 필요한 칸 수
            int local_Q1_before = Q1 % num_qubits;
            int dist_to_top_before = local_Q1_before / n;

            // 스왑 후: Q2 위치로 이동했을 때 최상단 Row까지 가기 위해 필요한 칸 수
            int local_Q2_after = Q2 % num_qubits;
            int dist_to_top_after = local_Q2_after / n;

            double current_cost = (dist_to_top_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost = (dist_to_top_after * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;


            final_cost += (current_cost - future_cost);
        }

        if (q2 == control || q2 == target) {

            int local_Q2_before = Q2 % num_qubits;
            int dist_to_top_before = local_Q2_before / n;

            int local_Q1_after = Q1 % num_qubits;
            int dist_to_top_after = local_Q1_after / n;

            double current_cost = (dist_to_top_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost = (dist_to_top_after * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;

            final_cost += (current_cost - future_cost);
        }
    }

    else //// Intra의 경우
    {
        // Cost: SWAP 시 두 대상 큐비트간의 거리

        final_cost = cal_MCPE(SWAP_pair, dgraph);
    }
    return final_cost;
}
/////////////////////////////////////////////////////////////////


void Qcircuit::QMapper::update_front_n_act_list(list<int>& front_list, list<int>& act_list, vector<bool>& frozen)
{
    for (int q = 0; q < nqubits; q++)
    {
        list<int>& Dlist_line = Dlist[q];
        if (!frozen[q])
        {
            if (Dlist_line.empty()) continue;
            int gateid = Dlist_line.front();

            if (find(front_list.begin(), front_list.end(), gateid) != front_list.end())
            {
                front_list.erase(remove(front_list.begin(), front_list.end(), gateid));
                act_list.push_back(gateid);
            }
            else
                front_list.push_back(gateid);
            frozen[q] = true;
        }
    }
}

void Qcircuit::QMapper::find_singlequbit_list(int gateid, list<int>& singlequbit_list, Circuit& dgraph)
{
    singlequbit_list.clear();
    int control = dgraph.nodeset[gateid].control;
    int target = dgraph.nodeset[gateid].target;

    int erase_control = 0;
    int erase_target = 0;

    //for control //
    for (auto& id : Dlist_all[control])
    {
        if (id == gateid) {
            erase_control++;
            break;
        }
        else {
            singlequbit_list.push_back(id);
            erase_control++;
        }
    }
    for (int i = 0; i < erase_control; i++) Dlist_all[control].pop_front();

    //for target //
    for (auto& id : Dlist_all[target])
    {
        if (id == gateid) {
            erase_target++;
            break;
        }
        else {
            singlequbit_list.push_back(id);
            erase_target++;
        }
    }
    for (int i = 0; i < erase_target; i++) Dlist_all[target].pop_front();

    singlequbit_list.sort();
}

void Qcircuit::QMapper::update_act_dist2_list(list<int>& act_dist2_list, list<int>& act_list, Circuit& dgraph)
{
    int idx1, idx2;
    act_dist2_list.clear();
    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target = dgraph.nodeset[gateid].target;
        idx1 = extract_qpu_idx(control);
        idx2 = extract_qpu_idx(target);

        bool is_inter_gate = (idx1 != idx2);
        if (is_inter_gate) continue; // 코어 간 연산인 경우 칩 경계를 넘는 스왑/브릿지는 차단

        int Q_control = layout_L[idx1][control];
        int Q_target = layout_L[idx2][target];
        if (multi_qpu_graph.dist[Q_control][Q_target] == 2)
            act_dist2_list.push_back(gateid);
    }
}

bool Qcircuit::QMapper::check_direct_act_list(list<int>& act_list, list<int>& singlequbit_list, vector<bool>& frozen, Circuit& dgraph)
{
    bool complete_act_list = false;
    vector<int> act_list_erase;

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target = dgraph.nodeset[gateid].target;

        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);

        int Q_control = layout_L[idx1][control];
        int Q_target = layout_L[idx2][target];

        bool can_execute = false;

        //////////////////////// bool 변수 업데이트 /////////////////////
        if (idx1 == idx2 && multi_qpu_graph.dist[Q_control][Q_target] == 1) {
            can_execute = true; // 같은 코어 내부 연산
        }

        else if (idx1 != idx2) {
            // 각자의 코어 내에서 로컬 위치가 최상단인지 여부
            bool control_on_top = (Q_control % num_qubits) < n;
            bool target_on_top = (Q_target % num_qubits) < n;

            if (control_on_top && target_on_top) {
                can_execute = true;
                // 여기서 나중에 EPR/Teleportation 게이트 등을 추가로 삽입할 수 있음
            }
        }
        /////////////////////////////////////////////////////////////////

        ///// 여기서 FinalCircuit에 추가함
        if (can_execute)
        {
            act_list_erase.push_back(gateid);
            Dlist[control].pop_front();
            Dlist[target].pop_front();
            frozen[control] = false;
            frozen[target] = false;

            FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
            complete_act_list = true;
        }
    }

    for (auto gateid : act_list_erase)
        act_list.erase(remove(act_list.begin(), act_list.end(), gateid));

    return complete_act_list;
}
///////////////////////////////////////////////////////////////////////필히 수정 요함, 인자에 플래그 넣어서 코어 안에 있는 연산만///////////////////////////////////////////////////////////////////////////////////////
void Qcircuit::QMapper::generate_candi_list(list<int>& act_list, vector< pair<pair<int, int>, int> >& candi_list, Circuit& dgraph, bool is_inter_gate)
{
    int idx1, idx2;

    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control; // 타깃 큐비트의 인덱스
        int target = dgraph.nodeset[gateid].target; // 제어 큐비트의 인덱스
        idx1 = extract_qpu_idx(control);
        idx2 = extract_qpu_idx(target);

        int Q_control = layout_L[idx1][control]; // 전체 하드웨어에서의 그 큐비트의 인덱스 말하는듯...?
        int Q_target = layout_L[idx2][target];

        bool is_inter_gate = (idx1 != idx2);

        for (int i = 0; i < multi_qpu_graph.node_size; i++)
        {   // 제어 큐비트에 대한 스왑 후보 탐색
            if (multi_qpu_graph.dist[Q_control][i] == 1)
            {
                if (extract_qpu_idx(i) != idx1) continue; // 코어 밖의 노드와는 스왑 불가능 하도록

                if (!is_inter_gate && cal_SWAP_effect(control, target, i, Q_control) <= 0) continue;

                (i < Q_control) ? candi_list.push_back(make_pair(make_pair(i, Q_control), gateid)) :
                    candi_list.push_back(make_pair(make_pair(Q_control, i), gateid));
            }
            // 타깃 큐비트에 대한 스왑 후보 탐색
            if (multi_qpu_graph.dist[Q_target][i] == 1)
            {
                if (extract_qpu_idx(i) != idx2) continue; // 코어 밖의 노드와는 스왑 불가능 하도록

                if (!is_inter_gate && cal_SWAP_effect(control, target, i, Q_target) <= 0) continue;
                (i < Q_target) ? candi_list.push_back(make_pair(make_pair(i, Q_target), gateid)) :
                    candi_list.push_back(make_pair(make_pair(Q_target, i), gateid));
            }
        }
        /////////////////////////////////////////////////////////////////////////
    }
}






void Qcircuit::QMapper::main_mapping(Circuit& dgraph, bool BRIDGE_MODE, int n)
{
    cout << "main_mapping\n";

    // (0) Make Dlist 
    make_Dlist(dgraph);
    make_Dlist_all(dgraph); // Dlist_all_mode이 0이 아닌 경우에만 동작

    // cost 관련 인자 선언
    add_2q_num = 0;
    fidelity = 0;

    ////////////////////////////////새롭게 추가한 COST 처리 관련 인자들/////////////////////////
    bool is_inter_gate = true; // 비용 함수 계산시  inter/intra 구분
    int loop_end = 0; // 전체 매인루프 종료 조건 인자
    //////////////////////////////////////////////////////////////////////////////////////////

    // (1) Circuit mapping 
    list<int> fron_list;
    list<int> act_list;
    list<int> act_dist2_list;
    vector<bool> frozen(nqubits, 0);
    list<int> singlequbit_list;

    FinalCircuit.nodeset.clear();
    node_id = dgraph.nodeset.size();
    add_cnot_num = 0;
    add_swap_num = 0;
    add_bridge_num = 0;

    bool in_process = true;

    vector< pair< pair<int, int>, pair<int, double> > > MCPE_flag;

    // ================= MAIN LOOP =================
    do {


        bool complete_act_list = true;

        do {
            update_front_n_act_list(fron_list, act_list, frozen);

            vector<int> buffer_usage(num_qpu, 0);
            vector<int> gates_to_defer;

            for (auto& gateid : act_list) {
                int control = dgraph.nodeset[gateid].control;
                int target = dgraph.nodeset[gateid].target;

                int core_c = extract_qpu_idx(layout_L[0][control]);
                int core_t = extract_qpu_idx(layout_L[0][target]);

                if (core_c != core_t) {
                    if (buffer_usage[core_c] < n && buffer_usage[core_t] < n) {
                        buffer_usage[core_c]++;
                        buffer_usage[core_t]++;
                    }
                    else {
                        gates_to_defer.push_back(gateid);
                    }
                }
            }

            for (auto& gateid : gates_to_defer) {
                act_list.remove(gateid);
                fron_list.push_back(gateid);

                int control = dgraph.nodeset[gateid].control;
                int target = dgraph.nodeset[gateid].target;
                frozen[control] = true;
                frozen[target] = true;
            }

            complete_act_list = check_direct_act_list(act_list, singlequbit_list, frozen, dgraph);
            act_list.sort();

        } while (complete_act_list);

        ////////////////////////////////두번째 do-while///////////////////////////////  
        do {

            vector< pair<pair<int, int>, int> > candi_list;
            pair<int, int> SWAP;
            int gateid;
            double max_cost;

            generate_candi_list(act_list, candi_list, dgraph, is_inter_gate);

            if (BRIDGE_MODE) {
                update_act_dist2_list(act_dist2_list, act_list, dgraph);
            }

            vector< pair< pair<int, int>, pair<int, double> > > MCPE_test;

            for (auto kv : candi_list) {
                pair<int, int> SWAP_pair = kv.first;
                int gateid = kv.second;

                int control = dgraph.nodeset[gateid].control;
                int target = dgraph.nodeset[gateid].target;

                bool is_inter_gate = (extract_qpu_idx(control) != extract_qpu_idx(target));
                double cost = mapping_machine(is_inter_gate, SWAP_pair, dgraph, gateid);

                MCPE_test.push_back(make_pair(SWAP_pair, make_pair(gateid, cost)));
            }

            find_max_cost(SWAP, gateid, max_cost, MCPE_test, act_list, MCPE_flag);

            do {

                if (!MCPE_test.empty()) {

                    int Q1 = SWAP.first;
                    int Q2 = SWAP.second;
                    int qpu_idx1 = Q1 / num_qubits;
                    int qpu_idx2 = Q2 / num_qubits;

                    int q1 = qubit_Q[qpu_idx1][Q1];
                    int q2 = qubit_Q[qpu_idx2][Q2];

                    int core1 = extract_qpu_idx(q1);
                    int core2 = extract_qpu_idx(q2);

                    layout_L[core1][q1] = Q2;
                    layout_L[core2][q2] = Q1;

                    qubit_Q[qpu_idx1][Q1] = q2;
                    qubit_Q[qpu_idx2][Q2] = q1;

                    add_swap(Q1, Q2, FinalCircuit);
                }

            } while (!act_list.empty());

            loop_end = 0;
            for (int q = 0; q < nqubits; q++) {
                if (Dlist[q].empty()) loop_end++;
            }
            if (loop_end == nqubits) in_process = false;

        } while (loop_end != nqubits); // ← 두번째 do 종료

    } while (in_process); // ← ★ 바깥 do 대응 추가
}
// ==================================================