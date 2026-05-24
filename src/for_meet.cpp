//for Meet
//main_mapping.cpp
#include "circuit.h"

using namespace std;
using namespace Qcircuit;

#define Dlist_all_mode 0 

void Qcircuit::QMapper::main_mapping(Circuit& dgraph){
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
    // 메인 매핑 진행하면서 게이트들을 분류할 리스트들
    list<int> fron_list;
    list<int> act_list;
    list<int> act_dist2_list;
    vector<bool> frozen(nqubits, 0);
    list<int> singlequbit_list;
    
    //initialize for post processing
    FinalCircuit.nodeset.clear(); 
    node_id = dgraph.nodeset.size();
    add_cnot_num = 0;
    add_swap_num = 0;
    add_bridge_num = 0;

    vector< pair< pair<int, int>, pair<int, double> > > MCPE_flag;

    // ================= MAIN LOOP =================
    do{

        bool complete_act_list = true;

        //////////////////////FSQM에서 그대로 가져오는 부분 + 수정////////////////////////////
        // DQC 환경에 맞게 update_front_n_act_list, check_direct_act_list 함수 수정해야함(커플링 그래프 도는 부분)
        do{
            // (1-1) Update front and act list
            update_front_n_act_list(fron_list, act_list, frozen);
            
            // inter 연산 시 하드웨어 제약 구현 부분
            vector<int> buffer_usage(num_qpu, 0); // 각 코어별로 현재 버퍼 카운트

            // 하드웨어 제약 구현 부분
            vector<int> gates_to_defer; // 버퍼 용량 초과로 이번 사이클에서 실행을 미룰 게이트 카운트

            for(auto& gateid : act_list) {
                int control = dgraph.nodeset[gateid].control;
                int target  = dgraph.nodeset[gateid].target;

                int core_c = extract_qpu_idx(layout_L[0][control]); // 제어 큐비트의 현재 코어
                int core_t = extract_qpu_idx(layout_L[0][target]);  // 타깃 큐비트의 현재 코어

                // Inter 연산(코어 간 통신)인 경우에만 버퍼를 차지함
                if(core_c != core_t) {
                    // 해당 코어들의 버퍼가 여유가 있는지 확인 (각각 n개 미만으로 사용 중이어야 함)
                    if(buffer_usage[core_c] < n && buffer_usage[core_t] < n) {
                        buffer_usage[core_c]++;
                        buffer_usage[core_t]++;
                    } else {
                        // 어느 한쪽이라도 버퍼 용량이 꽉 찼다면 이 게이트는 지금 실행 불가
                        gates_to_defer.push_back(gateid);
                    }
                }
            }

            // 미뤄진 게이트들 act_list -> fron_list + frozen
            for(auto& gateid : gates_to_defer) {
                act_list.erase(remove(act_list.begin(), act_list.end(), gateid));
                fron_list.push_back(gateid);
                
                int control = dgraph.nodeset[gateid].control;
                int target  = dgraph.nodeset[gateid].target;
                frozen[control] = true;
                frozen[target] = true;
            }

            // (1-2) Check direct act list
            complete_act_list = check_direct_act_list(act_list, singlequbit_list, frozen, dgraph);

            // (1-3) Sort act list
            act_list.sort();

        }while(complete_act_list);
        /////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////두번째 do-while///////////////////////////////  
        
    if(!act_lsit.empty()){

        // #1 후보 탐색(inter/intra 구분x)
        // Swap
        vector< pair<pair<int, int>, int> > candi_list;
        generate_candi_list(act_list, candi_list, dgraph); // 플래그 삽입

        // Bridge
        if(BRIDGE_MODE){
        update_act_dist2_list(act_dist2_list, act_list, dgraph); // 플래그 삽입
        }

        // #2 통합 Cost 계산
        vector< pair< pair<int, int>, pair<int, double> > > MCPE_test;
    
        for(auto kv : candi_list){
        pair<int, int> SWAP_pair = kv.first; // 가상의 SWAP 후보
        int gateid = kv.second;              // 이 SWAP이 해결하려는 타겟 게이트
        
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        
        // #3 gate 성격 판별
        bool is_inter_gate = (extract_qpu_idx(control) != extract_qpu_idx(target));

        // #4 mapping_machine으로 점수매기기 
        double cost = mapping_machine(is_inter_gate, SWAP_pair, dgraph, gateid); 
        
        MCPE_test.push_back(make_pair(SWAP_pair, make_pair(gateid, cost)));
        }

        int best_gateid = -1;
        double max_cost = -1e9; // 가장 낮은 값으로 초기화

        // 이 함수가 inter/intra 모두에게 다 돌아갈지 몰라서 코어 간 연산 비용 정의하고 수정해야할듯
        ////
        find_max_cost(SWAP, gateid, max_cost, MCPE_test, act_list, MCPE_flag);
        ////

         // #5 1등 후보 적용하고 레이아웃을 업데이트하는 기존 로직
        if (!MCPE_test.empty()){

            int Q1 = SWAP.first;
            int Q2 = SWAP.second;

            int q1 = qubit_Q[Q1];
            int q2 = qubit_Q[Q2];

            int core1 = extract_qpu_idx(q1) 
            int core2 = extract_qpu_idx(q2) 

            layout_L[core1][q1] = Q2;
            layout_L[core2][q2] = Q1;

            qubit_Q[Q1] = q2;
            qubit_Q[Q2] = q1;

            add_swap(Q1, Q2, FinalCircuit);
        }
    }
    /////////////////////////////////////////////////////////////////////////////

        // 메인 루프 종료 조건 갱신 (위치 고려 할 것)
        loop_end = 0; 
        for(int q=0; q<nqubits; q++){
            if(Dlist[q].empty()) loop_end++;
        }

    }while(loop_end != nqubits); /// 컴파일 되는지 보고 조건 수정해야함
}
// ==================================================
    
//////////////////////////// 추가한 함수 ////////////////////////
double Qcircuit::QMapper::mapping_machine(bool is_inter_gate, const pair<int, int> SWAP_pair, Circuit& dgraph, int gateid){
    
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;

    int q1 = qubit_Q[Q1];
    int q2 = qubit_Q[Q2];

    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    double final_cost = 0.0;
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    if (is_inter_gate) 
    {
        double WEIGHT_SWAP = 3.0; // 한 칸 전진할 때 얻는 점수(가중치)

        // 1. q1이 타깃/제어 큐비트일 때의 점수 계산
        if (q1 == control || q1 == target) {
            int dist_before = (Q1 % num_qubits) / n;
            int dist_after  = (Q2 % num_qubits) / n;
            // 상쇄되던 무의미한 버퍼 변수 제거. 순수 이동 점수만 계산
            final_cost += (dist_before - dist_after) * WEIGHT_SWAP; 
        }

        // 2. q2가 타깃/제어 큐비트일 때의 점수 계산
        if (q2 == control || q2 == target) {
            int dist_before = (Q2 % num_qubits) / n;
            int dist_after  = (Q1 % num_qubits) / n;
            final_cost += (dist_before - dist_after) * WEIGHT_SWAP;
        }

        // 3. 무의미한 횡이동(가로 SWAP) 페널티 (포인트 B 해결)
        // 위로 올라가지도, 내려가지도 않는 가로 이동이라면 미세한 감점(-0.1)을 주어 
        // 굳이 불필요한 SWAP이 선택되지 않도록 방어합니다.
        if (final_cost == 0.0) {
            final_cost -= 0.1; 
        }
    }
    else 
    {
        // Intra 연산의 경우: 기존 MCPE(거리 기반) 로직 그대로 수행
        final_cost = cal_MCPE(SWAP_pair, dgraph);
    }
    
    return final_cost;
}
/////////////////////////////////////////////////////////////////


void Qcircuit::QMapper::update_front_n_act_list(list<int>& front_list, list<int>& act_list, vector<bool>& frozen)
{
    for(int q = 0; q < nqubits; q++)
    {
        list<int>& Dlist_line = Dlist[q];
        if(!frozen[q])
        {
            if(Dlist_line.empty()) continue;
            int gateid = Dlist_line.front();
            
            if( find(front_list.begin(), front_list.end(), gateid) != front_list.end() )
            {
                front_list.erase( remove(front_list.begin(), front_list.end(), gateid));
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
    int target  = dgraph.nodeset[gateid].target;

    int erase_control = 0;
    int erase_target = 0;
    
    //for control //
    for(auto& id : Dlist_all[control])
    {
        if(id == gateid) {
            erase_control++;
            break;
        }
        else {
            singlequbit_list.push_back(id);
            erase_control++;
        }
    }
    for(int i=0; i<erase_control; i++) Dlist_all[control].pop_front();

    //for target //
    for(auto& id : Dlist_all[target])
    {
        if(id == gateid) {
            erase_target++;
            break;
        }
        else {
            singlequbit_list.push_back(id);
            erase_target++;
        }
    }
    for(int i=0; i<erase_target; i++) Dlist_all[target].pop_front();
    
    singlequbit_list.sort();
}

void Qcircuit::QMapper::update_act_dist2_list(list<int>& act_dist2_list, list<int>& act_list, Circuit& dgraph)
{
    int idx1,idx2;
    act_dist2_list.clear();
    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        idx1=extract_qpu_idx(control);
        idx2=extract_qpu_idx(target);

        bool is_inter_gate = (idx1 != idx2);
        if(is_inter_gate) continue; // 코어 간 연산인 경우 칩 경계를 넘는 스왑/브릿지는 차단

        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];
        if(multi_qpu_graph.dist[Q_control][Q_target] == 2)
            act_dist2_list.push_back(gateid);
    }
}

bool Qcircuit::QMapper::check_direct_act_list(list<int>& act_list, list<int>& singlequbit_list, vector<bool>& frozen, Circuit& dgraph)
{
    bool complete_act_list = false;
    vector<int> act_list_erase;

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    
    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;

        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);
        
        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];

        bool can_execute = false;

        //////////////////////// bool 변수 업데이트 /////////////////////
        if (idx1 == idx2 && multi_qpu_graph.dist[Q_control][Q_target] == 1) {
            can_execute = true; // 같은 코어 내부 연산
        }

        else if (idx1 != idx2) {
            // 각자의 코어 내에서 로컬 위치가 최상단인지 여부
            bool control_on_top = (Q_control % num_qubits) < n;
            bool target_on_top  = (Q_target % num_qubits) < n;

            if (control_on_top && target_on_top) {
                can_execute = true;
                // 여기서 나중에 EPR/Teleportation 게이트 등을 추가로 삽입할 수 있음
            }
        }
        /////////////////////////////////////////////////////////////////
        
        ///// 여기서 FinalCircuit에 추가함
        if(can_execute) 
        {
            act_list_erase.push_back(gateid);
            Dlist[control].pop_front();
            Dlist[target ].pop_front();
            frozen[control] = false;
            frozen[target ] = false;

            FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
            complete_act_list = true;
        }
    }
    
    for(auto gateid : act_list_erase)
        act_list.erase(remove(act_list.begin(), act_list.end(), gateid) );

    return complete_act_list;
}
///////////////////////////////////////////////////////////////////////필히 수정 요함, 인자에 플래그 넣어서 코어 안에 있는 연산만///////////////////////////////////////////////////////////////////////////////////////
void Qcircuit::QMapper::generate_candi_list(list<int>& act_list, vector< pair<pair<int, int>, int> >& candi_list, Circuit& dgraph, bool is_inter_gate)
{
    int idx1, idx2;

    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control; // 타깃 큐비트의 인덱스
        int target  = dgraph.nodeset[gateid].target; // 제어 큐비트의 인덱스
        idx1 = extract_qpu_idx(control);
        idx2 = extract_qpu_idx(target);
        
        int Q_control = layout_L[idx1][control]; // 전체 하드웨어에서의 그 큐비트의 인덱스 말하는듯...?
        int Q_target  = layout_L[idx2][target];

        bool is_inter_gate = (idx1 != idx2);
        
        for(int i = 0; i < multi_qpu_graph.node_size; i++)
        {   // 제어 큐비트에 대한 스왑 후보 탐색
            if(multi_qpu_graph.dist[Q_control][i] == 1)
            {
                if (extract_qpu_idx(i) != idx1) continue; // 코어 밖의 노드와는 스왑 불가능 하도록
                
                if(!is_inter_gate && cal_SWAP_effect(control, target, i, Q_control) <= 0) continue;

                (i < Q_control) ? candi_list.push_back(make_pair(make_pair(i, Q_control), gateid)) : 
                                  candi_list.push_back(make_pair(make_pair(Q_control, i), gateid));
            }
            // 타깃 큐비트에 대한 스왑 후보 탐색
            if(multi_qpu_graph.dist[Q_target][i] == 1)
            {
                if (extract_qpu_idx(i) != idx2) continue; // 코어 밖의 노드와는 스왑 불가능 하도록

                if(!is_inter_gate && cal_SWAP_effect(control, target, i, Q_target ) <= 0) continue;
                (i < Q_target)  ? candi_list.push_back(make_pair(make_pair(i, Q_target),  gateid)) : 
                                  candi_list.push_back(make_pair(make_pair(Q_target, i),  gateid));
            }
        }
        /////////////////////////////////////////////////////////////////////////
    }
}