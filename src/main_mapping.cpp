#include "circuit.h"
#include <algorithm>
#include <iostream>
#include <cmath>

using namespace std;
using namespace Qcircuit;

#define EMPTY_NODE -1  // 초기 상태 및 이니셜 매핑에서 선택되지 않은 빈 칸 정의

// ==========================================================
// 1. main_mapping() - 메인 제어 루프 (전체 재설계)
// ==========================================================
void Qcircuit::QMapper::main_mapping(
    Circuit& dgraph,
    bool BRIDGE_MODE,
    int n_buffer)
{
    cout << "[START] main_mapping - 빈 칸 활용 및 고립 방지 통합 버전\n";

    // 데이터 초기화 및 의존성 리스트 구축
    make_Dlist(dgraph);
    make_Dlist_all(dgraph);

    add_2q_num = 0;
    fidelity = 0;

    list<int> fron_list;
    list<int> act_list;
    list<int> singlequbit_list;

    // 각 물리 큐비트의 고정(락) 상태 관리
    vector<bool> frozen(nqubits, false);

    FinalCircuit.nodeset.clear();

    bool in_process = true;
    unsigned long long int safety_counter = 0; 

    do {
        bool complete_act_list = true;

        // --------------------------------------------------
        // STEP 1: 즉시 실행 가능한 게이트 선처리 및 업데이트
        // --------------------------------------------------
        do {
            update_front_n_act_list(fron_list, act_list, frozen);

            complete_act_list = check_direct_act_list(
                                    act_list,
                                    singlequbit_list,
                                    frozen,
                                    dgraph);

            act_list.unique(); // 중복 노드 제거 (순서는 유지)

        } while (complete_act_list);

        // --------------------------------------------------
        // STEP 2: 액티브 리스트 최전방 게이트 집중 라우팅
        // --------------------------------------------------
        if (!act_list.empty()) 
        {
            // 가장 의존성이 높고 시급한 맨 앞의 게이트 추출
            int gateid = act_list.front(); 

            int control = dgraph.nodeset[gateid].control;
            int target  = dgraph.nodeset[gateid].target;

            bool is_inter = (extract_qpu_idx(control) != extract_qpu_idx(target));
            bool executed = false;

            // 라우팅 프로토콜 실행 (act_list 참조 전달로 완결성 확보)
            if (is_inter) 
            {
                executed = handle_inter_gate(gateid, dgraph, frozen, act_list);
            } 
            else 
            {
                executed = handle_intra_gate(gateid, dgraph, frozen, act_list);
            }

            // 이번 루프에서 게이트가 실행되지 못했다면 다음 루프 탐색을 위해 락 해제
            if (!executed) 
            {
                frozen[control] = false;
                frozen[target]  = false;
            }
        }
        else
        {
            // 처리할 액티브 게이트가 일시적으로 없다면 전체 락을 풀고 대기
            fill(frozen.begin(), frozen.end(), false);
        }

        // --------------------------------------------------
        // STEP 3: 매핑 종료 조건 검사
        // --------------------------------------------------
        int loop_end = 0;
        for(int q = 0; q < nqubits; q++)
        {
            if(Dlist[q].empty())
                loop_end++;
        }

        if(loop_end == nqubits)
        {
            in_process = false;
        }

        // 🚨 최후의 대드락 보루 안전장치
        if (++safety_counter > 1000000) {
            add_2q_num=INT_MAX;
            break;
        }

    } while(in_process);

    cout << "[END] main_mapping 완료. 추가된 SWAP 개수: " << add_2q_num << "\n";
}

// ==========================================================
// 2. handle_inter_gate() - 글로벌 코어 간 라우팅
// ==========================================================
bool Qcircuit::QMapper::handle_inter_gate(
    int gateid,
    Circuit& dgraph,
    vector<bool>& frozen,
    list<int>& act_list)
{
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    int idx1 = extract_qpu_idx(control);
    int idx2 = extract_qpu_idx(target);

    int Qc = layout_L[idx1][control];
    int Qt = layout_L[idx2][target];

    int n = (int)std::sqrt(num_qubits);

    // 상단 통신 버퍼 진입 여부 체크
    bool c_top = ((Qc % num_qubits) < n);
    bool t_top = ((Qt % num_qubits) < n);

    // [성공] 둘 다 상단 버퍼에 도달했다면 인터 게이트 실행
    if(c_top && t_top)
    {
        FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);

        Dlist[control].pop_front();
        Dlist[target].pop_front();

        frozen[control] = false;
        frozen[target]  = false;

        act_list.remove(gateid); 
        return true;
    }

    // [진전] 한쪽만 안착한 경우, 안착한 쪽은 Hold하고 반대쪽만 전진 유도
    if(!c_top && t_top)
    {
        move_qubit_up(Qc);
        return false;
    }

    if(c_top && !t_top)
    {
        move_qubit_up(Qt);
        return false;
    }

    if(!c_top && !t_top)
    {
        move_qubit_up(Qc); // 순차적 라우팅을 위해 순서대로 전진
        return false;
    }

    return false;
}


bool Qcircuit::QMapper::handle_intra_gate(
    int gateid,
    Circuit& dgraph,
    vector<bool>& frozen,
    list<int>& act_list)
{
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    int idx1 = extract_qpu_idx(control);
    int idx2 = extract_qpu_idx(target);

    int Q_control = layout_L[idx1][control];
    int Q_target  = layout_L[idx2][target];

    // [성공] 이미 인접해 있다면 즉시 게이트 소화
    if(multi_qpu_graph.dist[Q_control][Q_target] == 1)
    {
        FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);

        Dlist[control].pop_front();
        Dlist[target].pop_front();

        frozen[control] = false;
        frozen[target]  = false;

        act_list.remove(gateid); 
        return true;
    }

    vector<pair<int,int>> candi;

    // 1차 탐색: 빈 칸 프리패스 및 점수가 유효한 스왑 후보군 수집 (force_mode = false)
    generate_intra_candidates(Q_control, Q_target, idx1, idx2, control, target, candi, false);

    // 2차 탐색: 주변에 유효한 진로가 없다면 코어 제약과 점수 제약을 전부 해제 (force_mode = true)
    if(candi.empty())
    {
        generate_intra_candidates(Q_control, Q_target, idx1, idx2, control, target, candi, true);
    }

    // 최후의 예외 처리 (그래프 상 이웃이 아예 없는 고립 상태 방지)
    if(candi.empty())
    {
        return false;
    }

    // 비용 함수 평가 후 최적의 스왑 적용
    double best_cost = -1e18;
    pair<int,int> best_swap = candi[0];

    for(auto& sw : candi)
    {
        double cost = cal_MCPE(sw, dgraph);
        if(cost > best_cost)
        {
            best_cost = cost;
            best_swap = sw;
        }
    }

    do_swap(best_swap.first, best_swap.second);
    return false;
}

// ==========================================================
// 4. generate_intra_candidates() - 후보군 생성기 (빈 칸 로직 탑재)
// ==========================================================
void Qcircuit::QMapper::generate_intra_candidates(
    int Q_control, int Q_target, int idx1, int idx2,
    int control, int target, vector<pair<int,int>>& candi, bool force_mode) 
{
    for(int i = 0; i < multi_qpu_graph.node_size; i++)
    {
        // --------------------------------------------------
        // Control 측 이웃 탐색
        // --------------------------------------------------
        if(multi_qpu_graph.dist[Q_control][i] == 1)
        {
            if(!force_mode && extract_qpu_idx(i) != idx1) continue; 

            int target_qpu = i / num_qubits;
            
            // 💡 [수정] map에 해당 물리 위치 i가 등록되어 있지 않다면(count == 0) 빈 칸입니다!
            if (qubit_Q[target_qpu].count(i) == 0) 
            {
                candi.push_back({ min(i, Q_control), max(i, Q_control) });
                continue; 
            }

            if(!force_mode && cal_SWAP_effect(control, target, i, Q_control) <= 0) 
                continue;

            candi.push_back({ min(i, Q_control), max(i, Q_control) });
        }

        // --------------------------------------------------
        // Target 측 이웃 탐색
        // --------------------------------------------------
        if(multi_qpu_graph.dist[Q_target][i] == 1)
        {
            if(!force_mode && extract_qpu_idx(i) != idx2) continue; 

            int target_qpu = i / num_qubits;

            // 💡 [수정] 동일하게 물리 위치 i가 map에 없으면 빈 칸으로 인정
            if (qubit_Q[target_qpu].count(i) == 0) 
            {
                candi.push_back({ min(i, Q_target), max(i, Q_target) });
                continue;
            }

            if(!force_mode && cal_SWAP_effect(control, target, i, Q_target) <= 0) 
                continue;

            candi.push_back({ min(i, Q_target), max(i, Q_target) });
        }
    }
}

// ==========================================================
// 5. do_swap() - 큐비트 위치 및 인덱스 맵 교환 (방어적 바인딩)
// ==========================================================
void Qcircuit::QMapper::do_swap(int Q1, int Q2)
{
    int qpu_idx1 = Q1 / num_qubits;
    int qpu_idx2 = Q2 / num_qubits;

    // 맵에 키가 있으면 해당 가상 큐비트 ID를 가져오고, 없으면 -1(임시 빈칸 표시) 처리
    int q1 = (qubit_Q[qpu_idx1].count(Q1) ? qubit_Q[qpu_idx1][Q1] : -1);
    int q2 = (qubit_Q[qpu_idx2].count(Q2) ? qubit_Q[qpu_idx2][Q2] : -1);

    // 1. layout_L (가상 -> 물리) 업데이트 및 기존 정보 삭제
    if (q1 != -1) {
        int core1 = extract_qpu_idx(q1);
        layout_L[core1][q1] = Q2;
    }
    if (q2 != -1) {
        int core2 = extract_qpu_idx(q2);
        layout_L[core2][q2] = Q1;
    }

    // 2. qubit_Q (물리 -> 가상) 교환 및 빈 칸 처리
    // Q1 자리에 q2 배치
    if (q2 != -1) {
        qubit_Q[qpu_idx1][Q1] = q2;
    } else {
        qubit_Q[qpu_idx1].erase(Q1); // q2가 빈 칸이었다면 Q1 위치는 이제 빈 칸이 됨
    }

    // Q2 자리에 q1 배치
    if (q1 != -1) {
        qubit_Q[qpu_idx2][Q2] = q1;
    } else {
        qubit_Q[qpu_idx2].erase(Q2); // q1이 빈 칸이었다면 Q2 위치는 이제 빈 칸이 됨
    }

    add_swap(Q1, Q2, FinalCircuit);
    add_2q_num++;
}
// ==========================================================
// 6. move_qubit_up() - 글로벌 버퍼 방향으로 전진
// ==========================================================
void Qcircuit::QMapper::move_qubit_up(int Q)
{
    int n = (int)std::sqrt(num_qubits);
    int qpu_idx = Q / num_qubits;
    int local = Q % num_qubits;
    int row = local / n;
    int col = local % n;

    if(row == 0) return; // 이미 최상단 버퍼열이면 스킵

    int upper_local = (row - 1) * n + col;
    int upper_global = qpu_idx * num_qubits + upper_local;

    do_swap(Q, upper_global);
}

// ==========================================================
// 7. update_front_n_act_list() - 전방 및 활성 의존성 리스트 갱신
// ==========================================================
void Qcircuit::QMapper::update_front_n_act_list(
    list<int>& front_list,
    list<int>& act_list,
    vector<bool>& frozen)
{
    for (int q = 0; q < nqubits; q++)
    {
        list<int>& Dlist_line = Dlist[q];

        if (frozen[q]) continue;
        if (Dlist_line.empty()) continue;

        int gateid = Dlist_line.front();
        auto front_it = find(front_list.begin(), front_list.end(), gateid);

        if (front_it != front_list.end())
        {
            front_list.remove(gateid);
            if(find(act_list.begin(), act_list.end(), gateid) == act_list.end())
            {
                act_list.push_back(gateid);
            }
        }
        else
        {
            front_list.push_back(gateid);
        }
        frozen[q] = true;
    }
}

// ==========================================================
// 8. check_direct_act_list() - 이미 붙어있는 인트라 게이트 즉시 실행
// ==========================================================
bool Qcircuit::QMapper::check_direct_act_list(
    list<int>& act_list,
    list<int>& singlequbit_list,
    vector<bool>& frozen,
    Circuit& dgraph)
{
    bool complete_act_list = false;
    vector<int> erase_list;

    for(auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;

        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);

        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];

        if(idx1 == idx2) 
        {
            if(multi_qpu_graph.dist[Q_control][Q_target] == 1)
            {
                erase_list.push_back(gateid);
                FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);

                Dlist[control].pop_front();
                Dlist[target].pop_front();

                frozen[control] = false;
                frozen[target]  = false;
                complete_act_list = true;
            }
        }
    }

    for(auto gateid : erase_list)
    {
        act_list.remove(gateid);
    }

    return complete_act_list;
}