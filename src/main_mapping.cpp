//main_mapping.cpp
#include "circuit.h"

using namespace std;
using namespace Qcircuit;

#define Dlist_all_mode 0 

void Qcircuit::QMapper::main_mapping(Circuit& dgraph){
    cout << "main_mapping\n";

    // (0) Make Dlist 
    make_Dlist(dgraph);
    make_Dlist_all(dgraph); // Dlist_all_mode이 0이 아닐때만 쓸모있을 듯

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
    
    //initialize for post processing
    FinalCircuit.nodeset.clear(); 
    node_id = dgraph.nodeset.size();
    add_cnot_num = 0;
    add_swap_num = 0;
    add_bridge_num = 0;

    vector< pair< pair<int, int>, pair<int, double> > > MCPE_flag;
    // int history_size = 2; -> 어떻게 활용해야할지

    // ================= MAIN LOOP =================
    do{

        bool complete_act_list = true;

        //////////////////////FSQM에서 그대로 가져오는 부분////////////////////////////
        // DQC 환경에 맞게 update_front_n_act_list, check_direct_act_list 함수 수정해야함(커플링 그래프 도는 부분)
        do{
            // (1-1) Update front and act list
            update_front_n_act_list(fron_list, act_list, frozen);

            // (1-2) Check direct act list
            complete_act_list = check_direct_act_list(act_list, singlequbit_list, frozen, dgraph);

            // (1-3) Sort act list
            act_list.sort();

        }while(complete_act_list);
        /////////////////////////////////////////////////////////////////////////////

        ////////////////////////////////두번째 do-while///////////////////////////////  
        do {

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

        // 이 함수가 inter/intra 모두에게 다 돌아갈지 몰라서 코어 간 연산 비용 정의하고 수정해야할듯
        find_max_cost(SWAP, gateid, max_cost, MCPE_test, act_list, MCPE_flag);

        // #5 1등 후보 적용하고 레이아웃을 업데이트하는 기존 로직
        //TO DO

        }while(!act_list.empty());
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
    
    // 가상 스왑 진행할 두 물큐의 인덱스
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;

    // 해당 물큐에 올라가 있는 논큐 인덱스
    int q1 = qubit_Q[Q1];
    int q2 = qubit_Q[Q2];

    // 타깃 게이트의 타깃, 제어 큐비트 인덱스
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    double final_cost = 0.0;

    int core_size = 9;  // 하드코딩 상태임
    int top_buffer_count = 3; // 시뮬마다 하드웨어 세팅에 맞게 바꿔주어야함.
    // 이거 관련 함수 만들기

    if (is_inter_gate)  //Inter의 경우
    {
        int min_dist_before = 9999;
        int min_dist_after = 9999;

        // Cost: SWAP 시 대상 큐비트와 버퍼 큐비트간의 거리

        // 그 논큐가 타코어와 연산해야하는 제어 또는 타깃인 경우
        if (q1 == control || q1 == target) {
            int core_idx = Q1 / core_size; // core_size 로 나누어 그 논큐의 코어 번호 획득
            
            // 유동적인 상단 버퍼 개수만큼 반복하며 어떤 버퍼로 가는 것이 가장빠른지 탐색
            for(int b = 0; b < top_buffer_count; b++) {
                int buf_node = (core_idx * core_size) + b; 
                //스왑 전 자리
                if(multi_qpu_graph.dist[Q1][buf_node] < min_dist_before) 
                    min_dist_before = multi_qpu_graph.dist[Q1][buf_node];
                //스왑 후 자리
                if(multi_qpu_graph.dist[Q2][buf_node] < min_dist_after) 
                    min_dist_after = multi_qpu_graph.dist[Q2][buf_node];
            }
            /* 최종 Cost 관련 멘트
            처음 거리 변수 설정 어케 할까.
            */

            final_cost += (min_dist_before - min_dist_after) * 10.0; // 10은 코어 간 연산
        }

        if (q2 == control || q2 == target) {
            int core_idx = Q2 / core_size; 

            for(int b = 0; b < top_buffer_count; b++) {
                int buf_node = (core_idx * core_size) + b;
                
                if(multi_qpu_graph.dist[Q2][buf_node] < min_dist_before) 
                    min_dist_before = multi_qpu_graph.dist[Q2][buf_node];
                    
                if(multi_qpu_graph.dist[Q1][buf_node] < min_dist_after) 
                    min_dist_after = multi_qpu_graph.dist[Q1][buf_node];
            }
            
            final_cost += (min_dist_before - min_dist_after) * 10.0;
        }
    }
    
    else //// Intra의 경우
    {
        // Cost: SWAP 시 두 대상 큐비트간의 거리

        final_cost = cal_MCPE(SWAP_pair, dgraph);

        // 물리적 노드 번호를 코어 크기로 나눈 나머지가 상단 가로 길이보다 작으면 버퍼 큐비트임
        bool is_Q1_buffer = ((Q1 % core_size) < top_buffer_count);
        bool is_Q2_buffer = ((Q2 % core_size) < top_buffer_count);

        // 내부 연산 큐비트가 상단 가장자리 버퍼 노드 쪽으로 침범하려 한다면 cost 조정
        if (is_Q1_buffer || is_Q2_buffer) {
            final_cost -= 500.0; // 이거 관련해서도,,,, 코스트 조정 잘 해야할듯......
            // 전반적인 cost가 너무 높거나 낮을 수도...? 근데 이게 상관있는 건지 없는 건지 모르겠음.
        }
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
    int idx1, idx2;
    vector<int> act_list_erase;
    
    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        idx1 = extract_qpu_idx(control);
        idx2 = extract_qpu_idx(target);
        
        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];
        
        if(multi_qpu_graph.dist[Q_control][Q_target] == 1) 
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
            
            //////////////////////////// 수정 ///////////////////////////////
            // 만약 코어 간 연산이라면 일반 CNOT 대신 EPR/Teleportation 게이트로 처리하도록
            if(idx1 != idx2){
            // TO DO
            }
            /////////////////////////////////////////////////////////////////

            FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
            complete_act_list = true;
        }
    }
    
    for(auto gateid : act_list_erase)
        act_list.erase( remove(act_list.begin(), act_list.end(), gateid) );

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