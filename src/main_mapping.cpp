// main_mapping.cpp
#include "circuit.h"
#include <climits>
#include <cmath>
#include <algorithm>
#include <iostream>

using namespace std;
using namespace Qcircuit;

#define Dlist_all_mode 0 

// BRIDGE_MODE가 매크로로 정의되어 있지 않을 경우를 대비한 안전장치
#ifndef BRIDGE_MODE
#define BRIDGE_MODE 0 
#endif

void Qcircuit::QMapper::main_mapping(Circuit& dgraph,bool BRIDGE_MODE){
    cout << "main_mapping\n";

    // (0) Make Dlist 
    make_Dlist(dgraph);
    make_Dlist_all(dgraph); // Dlist_all_mode이 0이 아닌 경우에만 동작

    // cost 관련 인자 선언
    add_2q_num = 0;
    fidelity = 0;

    ////////////////////////////////새롭게 추가한 COST 처리 관련 인자들/////////////////////////
    bool is_inter_gate = true; // 비용 함수 계산시  inter/intra 구분
    int loop_end = 0;          // 전체 매인루프 종료 조건 인자
    
    // 우우 // n과 num_qpu를 최상단에 선언하여 모든 루프 안에서 에러 없이 접근 가능하게 함
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    int num_qpu = multi_qpu_graph.node_size / num_qubits; 
    //////////////////////////////////////////////////////////////////////////////////////////

//     // (1) Circuit mapping 
//     list<int> fron_list;
//     list<int> act_list;
//     list<int> act_dist2_list;
//     vector<bool> frozen(nqubits, 0);
//     list<int> singlequbit_list;
    
    // initialize for post processing
    FinalCircuit.nodeset.clear(); 
    node_id = dgraph.nodeset.size();
    add_cnot_num = 0;
    add_swap_num = 0;
    add_bridge_num = 0;

    vector< pair< pair<int, int>, pair<int, double> > > MCPE_flag;

    // ================= MAIN LOOP =================
    do {
        bool complete_act_list = true;

        ////////////////////// FSQM에서 그대로 가져오는 부분 + 수정 ////////////////////////////
        do {
            // (1-1) Update front and act list
            update_front_n_act_list(fron_list, act_list, frozen);
            
            // inter 연산 시 하드웨어 제약 구현 부분
            vector<int> buffer_usage(num_qpu, 0); // 각 코어별로 현재 버퍼 카운트
            vector<int> gates_to_defer;           // 버퍼 용량 초과로 이번 사이클에서 실행을 미룰 게이트 카운트

            for(auto& gateid : act_list) {
                int control = dgraph.nodeset[gateid].control;
                int target  = dgraph.nodeset[gateid].target;

                // 우우 // 디버깅 논리 큐비트를 물리 큐비트로 변환(qubit_L) 후 코어 판별
                int core_c = extract_qpu_idx(qubit_L[control]); 
                int core_t = extract_qpu_idx(qubit_L[target]);  

                // Inter 연산(코어 간 통신)인 경우에만 버퍼를 차지함
                if(core_c != core_t) {
                    if(buffer_usage[core_c] < n && buffer_usage[core_t] < n) {
                        buffer_usage[core_c]++;
                        buffer_usage[core_t]++;
                    } else {
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

//             // (1-2) Check direct act list
//             complete_act_list = check_direct_act_list(act_list, singlequbit_list, frozen, dgraph);

//             // (1-3) Sort act list
//             act_list.sort();

        } while(complete_act_list);
        /////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////두번째 do-while///////////////////////////////  
        do {

            // #1 후보 탐색(inter/intra 구분x)
            vector< pair<pair<int, int>, int> > candi_list;
            generate_candi_list(act_list, candi_list, dgraph, is_inter_gate); // 플래그 인자 추가

            // Bridge
            if(BRIDGE_MODE){
                update_act_dist2_list(act_dist2_list, act_list, dgraph);
            }

            // #2 통합 Cost 계산
            vector< pair< pair<int, int>, pair<int, double> > > MCPE_test;
        
            for(auto kv : candi_list){
                pair<int, int> SWAP_pair = kv.first; 
                int loop_gateid = kv.second; // 변수명 충돌 방지
                
                int control = dgraph.nodeset[loop_gateid].control;
                int target  = dgraph.nodeset[loop_gateid].target;
                
                // #3 gate 성격 판별
                bool current_is_inter = (extract_qpu_idx(qubit_L[control]) != extract_qpu_idx(qubit_L[target]));

                // #4 mapping_machine으로 점수매기기 
                double cost = mapping_machine(current_is_inter, SWAP_pair, dgraph, loop_gateid); 
                
                MCPE_test.push_back(make_pair(SWAP_pair, make_pair(loop_gateid, cost)));
            }

            // 우우 /// find_max_cost 호출 전 파라미터 선언 필수
            pair<int, int> SWAP;
            int target_gateid; 
            double max_cost = 0.0;
            find_max_cost(SWAP, target_gateid, max_cost, MCPE_test, act_list, MCPE_flag);

            // #5 1등 후보 적용하고 레이아웃을 업데이트하는 기존 로직
            do {
                if (!MCPE_test.empty()){
                    int Q1 = SWAP.first;
                    int Q2 = SWAP.second;

                    int q1 = qubit_Q[Q1];
                    int q2 = qubit_Q[Q2];

                    // 우우 // 물리 인덱스 전달 및 세미콜론 추가
                    int core1 = extract_qpu_idx(Q1); 
                    int core2 = extract_qpu_idx(Q2); 

                    layout_L[core1][q1] = Q2;
                    layout_L[core2][q2] = Q1;

                    qubit_Q[Q1] = q2;
                    qubit_Q[Q2] = q1;
                    
                    // 물리적 위치 변동 시 qubit_L 업데이트 
                    qubit_L[q1] = Q2;
                    qubit_L[q2] = Q1;

                    add_swap(Q1, Q2, FinalCircuit);
                }
            } while(!act_list.empty());

            // 메인 루프 종료 조건 갱신 
            loop_end = 0; 
            for(int q=0; q<nqubits; q++){
                if(Dlist[q].empty()) loop_end++;
            }

        } while(loop_end != nqubits); 
    } while(loop_end != nqubits); 
}
// ==================================================
    
//////////////////////////// 추가한 함수 ////////////////////////
double Qcircuit::QMapper::mapping_machine(bool is_inter_gate, const pair<int, int> SWAP_pair, Circuit& dgraph, int gateid){
    
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;

    int q1 = qubit_Q[Q1];
    int q2 = qubit_Q[Q2];

//     int q1 = qubit_Q[idx1][Q1];
//     int q2 = qubit_Q[idx2][Q2];

//     int control = dgraph.nodeset[gateid].control;
//     int target  = dgraph.nodeset[gateid].target;

    double final_cost = 0.0;

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    double WEIGHT_SWAP = 3.0;           
    double WEIGHT_VIRTUAL_BUFFER = 10.0; 

    if (is_inter_gate) 
    {
        if (q1 == control || q1 == target) {
            int local_Q1_before = Q1 % num_qubits; 
            int dist_to_top_before = local_Q1_before / n; 
            
            int local_Q2_after = Q2 % num_qubits; 
            int dist_to_top_after = local_Q2_after / n; 

            double current_cost = (dist_to_top_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost  = (dist_to_top_after * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
    
            final_cost += (current_cost - future_cost);
        }

        if (q2 == control || q2 == target) {
            int local_Q2_before = Q2 % num_qubits;
            int dist_to_top_before = local_Q2_before / n;

            int local_Q1_after = Q1 % num_qubits;
            int dist_to_top_after = local_Q1_after / n;

            double current_cost = (dist_to_top_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost  = (dist_to_top_after * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;

            final_cost += (current_cost - future_cost);
        }
    }
    else 
    {
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
        
        idx1=extract_qpu_idx(qubit_L[control]);
        idx2=extract_qpu_idx(qubit_L[target]);

        bool is_inter_gate = (idx1 != idx2);
        if(is_inter_gate) continue; 

        int Q_control = qubit_L[control];
        int Q_target  = qubit_L[target];
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

        //우우// 물리 큐비트 기준으로 정확한 코어 인덱스 추출
        int idx1 = extract_qpu_idx(qubit_L[control]);
        int idx2 = extract_qpu_idx(qubit_L[target]);
        
        int Q_control = qubit_L[control];
        int Q_target  = qubit_L[target];

        bool can_execute = false;

        //////////////////////// bool 변수 업데이트 /////////////////////
        if (idx1 == idx2 && multi_qpu_graph.dist[Q_control][Q_target] == 1) {
            can_execute = true; // 같은 코어 내부 연산
        }
        else if (idx1 != idx2) {
            bool control_on_top = (Q_control % num_qubits) < n;
            bool target_on_top  = (Q_target % num_qubits) < n;

            if (control_on_top && target_on_top) {
                can_execute = true;
            }
        }
        /////////////////////////////////////////////////////////////////
        
        if(can_execute) 
        {
            act_list_erase.push_back(gateid);
            Dlist[control].pop_front();
            Dlist[target ].pop_front();
            frozen[control] = false;
            frozen[target ] = false;

            #if Dlist_all_mode
            find_singlequbit_list(gateid, singlequbit_list, dgraph);
            for(auto& id : singlequbit_list)
                FinalCircuit.nodeset.push_back(dgraph.nodeset[id]);
            #endif

            FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
            complete_act_list = true;
        }
    }
    
    for(auto gateid : act_list_erase)
        act_list.erase(remove(act_list.begin(), act_list.end(), gateid) );

    return complete_act_list;
}

void Qcircuit::QMapper::generate_candi_list(list<int>& act_list, vector< pair<pair<int, int>, int> >& candi_list, Circuit& dgraph, bool is_inter_gate)
{
    int idx1, idx2;

    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control; 
        int target  = dgraph.nodeset[gateid].target; 
        
        idx1 = extract_qpu_idx(qubit_L[control]);
        idx2 = extract_qpu_idx(qubit_L[target]);
        
        int Q_control = qubit_L[control]; 
        int Q_target  = qubit_L[target];

        bool current_is_inter = (idx1 != idx2);
        
        for(int i = 0; i < multi_qpu_graph.node_size; i++)
        {   
            if(multi_qpu_graph.dist[Q_control][i] == 1)
            {
                if (extract_qpu_idx(i) != idx1) continue; 
                
                if(!current_is_inter && cal_SWAP_effect(control, target, i, Q_control) <= 0) continue;

                (i < Q_control) ? candi_list.push_back(make_pair(make_pair(i, Q_control), gateid)) : 
                                  candi_list.push_back(make_pair(make_pair(Q_control, i), gateid));
            }
            
            if(multi_qpu_graph.dist[Q_target][i] == 1)
            {
                if (extract_qpu_idx(i) != idx2) continue; 

                if(!current_is_inter && cal_SWAP_effect(control, target, i, Q_target ) <= 0) continue;
                (i < Q_target)  ? candi_list.push_back(make_pair(make_pair(i, Q_target),  gateid)) : 
                                  candi_list.push_back(make_pair(make_pair(Q_target, i),  gateid));
            }
        }
    }
}