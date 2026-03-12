//main_mapping.cpp
#include "circuit.h"

using namespace std;
using namespace Qcircuit;

// < 표식 정리>
// Q : 어떤 기능인지 정확히 파악하지 못함
// R : 수정해야할 부분
// E : FSQM 그대로 가져갈 부분

#define Dlist_all_mode 0 // 1Q의 부분도 스케쥴링 할 것인지 결정하는 것 TEST를 위해선 0으로 놓고 하는게 좋을듯

void Qcircuit::QMapper::main_mapping(Circuit& dgraph){
    cout << "main_mapping\n";

    // (0) Make Dlist 
    make_Dlist(dgraph);
    make_Dlist_all(dgraph);

    add_2q_num = 0;
    fidelity = 0;

    ////////////////////////////////새롭게 추가한 COST 처리 관련 인자들/////////////////////////
     bool cost_flag = true; // 비용 함수 계산시  inter/intra 구분
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
        generate_candi_list(act_list, candi_list, dgraph);
        // Bridge
        if(BRIDGE_MODE){
        update_act_dist2_list(act_dist2_list, act_list, dgraph);
        }

        // #2 통합 Cost 계산
        vector< pair< pair<int, int>, pair<int, double> > > MCPE_test;
    
        for(auto kv : candi_list){
        pair<int, int> SWAP_pair = kv.first; // 가상의 SWAP 후보
        int gateid = kv.second;              // 이 SWAP이 해결하려는 타겟 게이트
        
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        
        // #3 gate 성격 판별
        bool cost_flag = (extract_qpu_idx(control) != extract_qpu_idx(target));

        // #4 mapping_machine으로 점수매기기 
        double cost = mapping_machine(cost_flag, SWAP_pair, dgraph, gateid); 
        
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
// Qmapper에 선언 추가해야함
void Qcircuit::QMapper::extract_qpu_idx_from_node(){
    // TO DO
    // 물리적 노드 번호 i를 넣었을 때 그 노드가 몇 번째 QPU 코어에 속해있는지 반환하도록
    // 기존 자료구조로 통합해야할듯
    // 무한루프 방지로 만듬
}
/////////////////////////////////////////////////////////////////
    
//////////////////////////// 추가한 함수 ////////////////////////
double Qcircuit::QMapper::mapping_machine(bool cost_flag, const pair<int, int> SWAP_pair, Circuit& dgraph, int gateid){
    
    // SWAP을 가상으로 진행해 볼 두 물리적 노드
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;

    // 해당 노드에 올라가 있는 논큐 인덱스
    int q1 = qubit_Q[Q1];
    int q2 = qubit_Q[Q2];

    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    double final_cost = 0.0;

    if (cost_flag)  //== Inter의 경우(코어간) ==
    {
        // 큐비트간 거리 저장 변수 선언
        int dist_before = 0;
        int dist_after = 0;

        // Cost: SWAP 시 대상 큐비트와 버퍼 큐비트간의 거리

        /* TO DO

        1. 해당 코어의 버퍼 큐비트 인덱스, 해당 큐비트의 인덱스 리턴

        2. q1과 버퍼 큐비트들 중 가장 짧은 거리, SWAP 후 Q1의 위치와 버퍼 큐비트들 중 가장 짧은 거리

        3. 반복문 돌면서 게이트 자료구조 돌면서 스왑전, 스왑후 거리 업데이트

        # 선택 요소  미래 게이트(Look-ahead) 평가 // 복잡한 요소
        q1, q2의 미래 스케줄을 보고, 미래에도 Inter 연산이 있다면 버퍼 근처에 머무는 것에 가중치 부여
        
        final_cost += ...;
        
        */
    }
    
    else ////== Intra의 경우(코어안) ==
    {
        // Cost: SWAP 시 두 대상 큐비트간의 거리
        
        /* TO DO

        1. 기존 FSQM의 cal_MCPE 로직을 활용하여 서로 거리가 줄어드는지 평가(Cal_MCPA 활용)
        
        2. 코어안 연산에서 버퍼 큐비트로 안가도록 가중치 설정

        final_cost += ...;
        */
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
        if(idx1 != idx2) continue;
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

void Qcircuit::QMapper::generate_candi_list(list<int>& act_list, vector< pair<pair<int, int>, int> >& candi_list, Circuit& dgraph)
{
    int idx1, idx2;
    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        idx1 = extract_qpu_idx(control);
        idx2 = extract_qpu_idx(target);
        
        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];
        
        //////////////////////////// 무한 루프 수정 ///////////////////////////////
        for(int i = 0; i < multi_qpu_graph.node_size; i++)
        {   // SWAP 대상 노드 i가 제어 큐비트와 같은 코어(idx1) 안에 있어야 조건 충족
            if(multi_qpu_graph.dist[Q_control][i] == 1 && extract_qpu_idx_from_node(i) == idx1)
            {
                // 코어가 다를 경우 cal_SWAP_effect는 상대방 코어의 경계선 노드와의 거리가 줄었는지 평가
                if(cal_SWAP_effect(control, target, i, Q_control) <= 0) continue;
                (i < Q_control) ? candi_list.push_back(make_pair(make_pair(i, Q_control), gateid)) : 
                                  candi_list.push_back(make_pair(make_pair(Q_control, i), gateid));
            }
            // SWAP 대상 노드 i가 대상 큐비트와 같은 코어(idx2) 안에 잇어야 조건 충족
            if(multi_qpu_graph.dist[Q_target][i] == 1 && extract_qpu_idx_from_node(i) == idx2)
            {
                if(cal_SWAP_effect(control, target, i, Q_target ) <= 0) continue;
                (i < Q_target)  ? candi_list.push_back(make_pair(make_pair(i, Q_target),  gateid)) : 
                                  candi_list.push_back(make_pair(make_pair(Q_target, i),  gateid));
            }
        }
        /////////////////////////////////////////////////////////////////////////
    }
}