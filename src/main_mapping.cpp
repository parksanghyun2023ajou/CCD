#include "circuit.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <climits>

using namespace std;
using namespace Qcircuit;

#define EMPTY_NODE -1

// =====================================================================
// 헬퍼: 빈 칸 안전한 논리큐비트 조회
// qubit_Q[qpu_idx]에 키가 없으면 -1 반환 (map의 [] 자동 삽입 방지)
// =====================================================================
static inline int safe_get_logical(
    vector<map<int,int>>& qubit_Q,
    int qpu_idx,
    int physical)
{
    auto it = qubit_Q[qpu_idx].find(physical);
    if (it == qubit_Q[qpu_idx].end()) return -1;
    return it->second;
}

// =====================================================================
// 헬퍼: 큐비트가 상단 row(버퍼 진입 가능 위치)에 있는지 확인
// =====================================================================
static inline bool is_on_top_row(int Q, int num_qubits, int n)
{
    int local = Q % num_qubits;
    return (local / n) == 0;  // row 0
}

// =====================================================================
// 헬퍼: 큐비트의 상단 row까지 거리 (몇 칸 위로 올라가야 하는지)
// =====================================================================
static inline int dist_to_top(int Q, int num_qubits, int n)
{
    int local = Q % num_qubits;
    return local / n;  // row 번호 = 위로 올라갈 칸 수
}


// =====================================================================
// 1. mapping_machine - 통합 비용 함수
//    네 원본 철학을 유지하되 빈 칸 안전 처리 추가
// =====================================================================
double Qcircuit::QMapper::mapping_machine(
    bool is_inter_gate,
    pair<int,int> SWAP_pair,
    Circuit& dgraph,
    int gateid)
{
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;
    int qpu_idx1 = Q1 / num_qubits;
    int qpu_idx2 = Q2 / num_qubits;

    // 빈 칸 안전 조회
    int q1 = safe_get_logical(qubit_Q, qpu_idx1, Q1);
    int q2 = safe_get_logical(qubit_Q, qpu_idx2, Q2);

    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    // 가중치
    const double WEIGHT_SWAP = 3.0;
    const double WEIGHT_VIRTUAL_BUFFER = 10.0;

    if (is_inter_gate)
    {
        // === Inter-gate: 상단 버퍼까지 거리 감소량 평가 ===
        double final_cost = 0.0;

        // q1이 inter-게이트의 핵심 큐비트인 경우
        if (q1 != -1 && (q1 == control || q1 == target)) {
            int dist_before = dist_to_top(Q1, num_qubits, n);
            int dist_after  = dist_to_top(Q2, num_qubits, n);

            double current_cost = (dist_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost  = (dist_after  * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;

            // 거리가 줄어드는 쪽이 양수 cost
            final_cost += (current_cost - future_cost);
        }

        // q2가 inter-게이트의 핵심 큐비트인 경우
        if (q2 != -1 && (q2 == control || q2 == target)) {
            int dist_before = dist_to_top(Q2, num_qubits, n);
            int dist_after  = dist_to_top(Q1, num_qubits, n);

            double current_cost = (dist_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost  = (dist_after  * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;

            final_cost += (current_cost - future_cost);
        }

        // 빈 칸과의 SWAP은 자유롭게 허용 (cost = 0이 아니라 약한 양수)
        // → 빈 칸을 활용한 우회 경로 가능
        if ((q1 == -1 && q2 != -1 && (q2 == control || q2 == target)) ||
            (q2 == -1 && q1 != -1 && (q1 == control || q1 == target)))
        {
            // 핵심 큐비트가 빈 칸으로 이동 → 자유 이동 보너스
            if (final_cost == 0.0) final_cost = 1.0;
        }

        return final_cost;
    }
    else
    {
        // === Intra-gate: MCPE 비용 사용 (네 원본 철학) ===
        // 단, 빈 칸이 끼면 MCPE가 무의미하므로 작은 양수 부여
        if (q1 == -1 || q2 == -1) {
            return 0.5;  // 빈 칸 SWAP 약한 허용
        }
        return cal_MCPE(SWAP_pair, dgraph);
    }
}


// =====================================================================
// 2. update_front_n_act_list - 전방/활성 리스트 갱신
//    (친구 코드 기반, 약간 개선)
// =====================================================================
void Qcircuit::QMapper::update_front_n_act_list(
    list<int>& front_list,
    list<int>& act_list,
    vector<bool>& frozen)
{
    for (int q = 0; q < (int)nqubits; q++)
    {
        if (frozen[q]) continue;
        if (Dlist[q].empty()) continue;

        int gateid = Dlist[q].front();

        auto front_it = find(front_list.begin(), front_list.end(), gateid);
        if (front_it != front_list.end())
        {
            // 두 큐비트 모두 같은 게이트를 front로 가지면 act로 승격
            front_list.erase(front_it);
            if (find(act_list.begin(), act_list.end(), gateid) == act_list.end())
                act_list.push_back(gateid);
        }
        else
        {
            front_list.push_back(gateid);
        }
        frozen[q] = true;
    }
}


// =====================================================================
// 3. check_direct_act_list - 즉시 실행 가능한 게이트 처리
//    (intra: 인접 / inter: 양쪽 모두 top row)
// =====================================================================
bool Qcircuit::QMapper::check_direct_act_list(
    list<int>& act_list,
    list<int>& singlequbit_list,
    vector<bool>& frozen,
    Circuit& dgraph)
{
    bool complete = false;
    vector<int> erase_list;

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;

        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);

        // layout_L 안전 조회
        if (layout_L[idx1].count(control) == 0) continue;
        if (layout_L[idx2].count(target)  == 0) continue;

        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];

        bool can_execute = false;

        if (idx1 == idx2) {
            // Intra: 그래프 상 거리 1이면 즉시 실행
            if (multi_qpu_graph.dist[Q_control][Q_target] == 1)
                can_execute = true;
        }
        else {
            // Inter: 양쪽 모두 상단 row면 실행 (가상 버퍼 통신)
            bool c_top = is_on_top_row(Q_control, num_qubits, n);
            bool t_top = is_on_top_row(Q_target,  num_qubits, n);
            if (c_top && t_top)
                can_execute = true;
        }

        if (can_execute)
        {
            erase_list.push_back(gateid);
            FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);

            Dlist[control].pop_front();
            Dlist[target].pop_front();
            frozen[control] = false;
            frozen[target]  = false;
            complete = true;
        }
    }

    for (auto gid : erase_list)
        act_list.remove(gid);

    return complete;
}


// =====================================================================
// 4. do_swap - 빈 칸 안전한 SWAP 적용
//    (친구 코드의 방어적 처리 유지 + 일관성 강화)
// =====================================================================
void Qcircuit::QMapper::do_swap(int Q1, int Q2)
{
    int qpu_idx1 = Q1 / num_qubits;
    int qpu_idx2 = Q2 / num_qubits;

    int q1 = safe_get_logical(qubit_Q, qpu_idx1, Q1);
    int q2 = safe_get_logical(qubit_Q, qpu_idx2, Q2);

    // 1) layout_L 업데이트 (논리 → 물리)
    if (q1 != -1) {
        int core1 = extract_qpu_idx(q1);
        if (core1 >= 0) layout_L[core1][q1] = Q2;
    }
    if (q2 != -1) {
        int core2 = extract_qpu_idx(q2);
        if (core2 >= 0) layout_L[core2][q2] = Q1;
    }

    // 2) qubit_Q 업데이트 (물리 → 논리), 빈 칸은 erase로 명확화
    if (q2 != -1)
        qubit_Q[qpu_idx1][Q1] = q2;
    else
        qubit_Q[qpu_idx1].erase(Q1);

    if (q1 != -1)
        qubit_Q[qpu_idx2][Q2] = q1;
    else
        qubit_Q[qpu_idx2].erase(Q2);

    // 3) 실제 SWAP 게이트 추가는 두 큐비트 중 하나라도 의미 있는 경우만
    //    (빈 칸끼리의 SWAP은 회로에 추가할 필요 없음)
    if (q1 != -1 || q2 != -1)
    {
        add_swap(Q1, Q2, FinalCircuit);
        add_2q_num++;
    }
}


// =====================================================================
// 5. generate_intra_candidates - 같은 코어 내 SWAP 후보 생성
// =====================================================================
void Qcircuit::QMapper::generate_intra_candidates(
    int Q_control, int Q_target,
    int idx1, int idx2,
    int control, int target,
    vector<pair<int,int>>& candi,
    bool force_mode)
{
    for (int i = 0; i < multi_qpu_graph.node_size; i++)
    {
        // --- Control 이웃 ---
        if (multi_qpu_graph.dist[Q_control][i] == 1)
        {
            if (!force_mode && extract_qpu_idx(i) != idx1 && (i / num_qubits) != idx1) {
                // 코어 밖 후보는 force_mode에서만 허용
                // (extract_qpu_idx는 logical 기준이므로 i가 빈 칸일 수 있어 추가 체크)
            }
            else
            {
                int target_qpu = i / num_qubits;
                bool is_empty = (qubit_Q[target_qpu].count(i) == 0);

                if (is_empty) {
                    // 빈 칸은 free SWAP (이동만 발생)
                    candi.push_back({min(i, Q_control), max(i, Q_control)});
                }
                else if (force_mode ||
                         cal_SWAP_effect(control, target, i, Q_control) > 0)
                {
                    candi.push_back({min(i, Q_control), max(i, Q_control)});
                }
            }
        }

        // --- Target 이웃 ---
        if (multi_qpu_graph.dist[Q_target][i] == 1)
        {
            if (!force_mode && extract_qpu_idx(i) != idx2 && (i / num_qubits) != idx2) {
                // 위와 동일
            }
            else
            {
                int target_qpu = i / num_qubits;
                bool is_empty = (qubit_Q[target_qpu].count(i) == 0);

                if (is_empty) {
                    candi.push_back({min(i, Q_target), max(i, Q_target)});
                }
                else if (force_mode ||
                         cal_SWAP_effect(control, target, i, Q_target) > 0)
                {
                    candi.push_back({min(i, Q_target), max(i, Q_target)});
                }
            }
        }
    }

    // 중복 제거
    sort(candi.begin(), candi.end());
    candi.erase(unique(candi.begin(), candi.end()), candi.end());
}


// =====================================================================
// 6. generate_inter_candidates - inter-게이트용 위로 가는 후보 생성
//    핵심: "위쪽 방향(상단 row 쪽)"으로 가는 SWAP만 후보로 채택
//    → 친구 코드의 move_qubit_up을 비용 함수와 결합한 형태
// =====================================================================
static void generate_inter_candidates(
    int Q_focus,
    int num_qubits,
    int n,
    Graph& multi_qpu_graph,
    vector<map<int,int>>& qubit_Q,
    vector<pair<int,int>>& candi)
{
    int qpu_idx = Q_focus / num_qubits;
    int local = Q_focus % num_qubits;
    int row = local / n;
    int col = local % n;

    if (row == 0) return;

    vector<int> neighbors;

    // 위쪽 (필수)
    int up = qpu_idx * num_qubits + (row - 1) * n + col;
    neighbors.push_back(up);

    // 대각 위-좌
    if (col > 0) {
        int up_left = qpu_idx * num_qubits + (row - 1) * n + (col - 1);
        if (multi_qpu_graph.dist[Q_focus][up_left] == 1)
            neighbors.push_back(up_left);
    }
    // 대각 위-우
    if (col < n - 1) {
        int up_right = qpu_idx * num_qubits + (row - 1) * n + (col + 1);
        if (multi_qpu_graph.dist[Q_focus][up_right] == 1)
            neighbors.push_back(up_right);
    }

    for (int pos : neighbors)
    {
        if (pos / num_qubits != qpu_idx) continue;
        if (multi_qpu_graph.dist[Q_focus][pos] != 1) continue;
        candi.push_back({min(Q_focus, pos), max(Q_focus, pos)});
    }
}


// =====================================================================
// 7. handle_intra_gate - 같은 코어 내 게이트 처리
//    네 mapping_machine을 비용으로 사용
// =====================================================================
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

    if (layout_L[idx1].count(control) == 0) return false;
    if (layout_L[idx2].count(target)  == 0) return false;

    int Q_control = layout_L[idx1][control];
    int Q_target  = layout_L[idx2][target];

    // 이미 인접 → 즉시 실행
    if (multi_qpu_graph.dist[Q_control][Q_target] == 1)
    {
        FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
        Dlist[control].pop_front();
        Dlist[target].pop_front();
        frozen[control] = false;
        frozen[target]  = false;
        act_list.remove(gateid);
        return true;
    }

    // 후보 생성 (1차: 점수 양수만)
    vector<pair<int,int>> candi;
    generate_intra_candidates(Q_control, Q_target, idx1, idx2,
                              control, target, candi, false);

    // 2차: 강제 모드 (점수 무시)
    if (candi.empty())
        generate_intra_candidates(Q_control, Q_target, idx1, idx2,
                                  control, target, candi, true);

    if (candi.empty()) return false;

    // 비용 평가 → 최적 SWAP 선택
    double best_cost = -1e18;
    pair<int,int> best_swap = candi[0];

    for (auto& sw : candi)
    {
        double cost = mapping_machine(false, sw, dgraph, gateid);
        if (cost > best_cost)
        {
            best_cost = cost;
            best_swap = sw;
        }
    }

    do_swap(best_swap.first, best_swap.second);
    return false;  // 게이트는 아직 실행 안 됨 (다음 루프에서 재평가)
}


// =====================================================================
// 8. handle_inter_gate - 코어 간 게이트 처리
//    네 mapping_machine 비용으로 "어느 방향이 더 좋은지" 평가
// =====================================================================
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

    if (layout_L[idx1].count(control) == 0) return false;
    if (layout_L[idx2].count(target)  == 0) return false;

    int Qc = layout_L[idx1][control];
    int Qt = layout_L[idx2][target];

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    bool c_top = is_on_top_row(Qc, num_qubits, n);
    bool t_top = is_on_top_row(Qt, num_qubits, n);

    // [완료] 둘 다 상단 → inter-gate 실행 (가상 버퍼 통신)
    if (c_top && t_top)
    {
        FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
        Dlist[control].pop_front();
        Dlist[target].pop_front();
        frozen[control] = false;
        frozen[target]  = false;
        act_list.remove(gateid);
        return true;
    }

    // 위로 올려야 할 큐비트들 식별
    vector<int> need_up;
    if (!c_top) need_up.push_back(Qc);
    if (!t_top) need_up.push_back(Qt);

    // 후보 생성: 각각의 큐비트가 위로 가는 SWAP들
    vector<pair<int,int>> candi;
    for (int Q : need_up)
        generate_inter_candidates(Q, num_qubits, n,
                                  multi_qpu_graph, qubit_Q, candi);

    // 중복 제거
    sort(candi.begin(), candi.end());
    candi.erase(unique(candi.begin(), candi.end()), candi.end());

    if (candi.empty()) return false;

    // 비용 평가 → 최적 SWAP 선택
    double best_cost = -1e18;
    pair<int,int> best_swap = candi[0];

    for (auto& sw : candi)
    {
        double cost = mapping_machine(true, sw, dgraph, gateid);
        if (cost > best_cost)
        {
            best_cost = cost;
            best_swap = sw;
        }
    }

    // best_cost가 음수면 (오히려 멀어지는 SWAP만 있는 경우) 진행 안 함
    // → 데드락 회피
    if (best_cost <= 0) return false;

    do_swap(best_swap.first, best_swap.second);
    return false;
}


// =====================================================================
// 9. main_mapping - 메인 루프 (통합 버전)
// =====================================================================
void Qcircuit::QMapper::main_mapping(
    Circuit& dgraph,
    bool BRIDGE_MODE,
    int n_buffer)
{
    cout << "[START] main_mapping (Integrated Version)\n";

    make_Dlist(dgraph);
    make_Dlist_all(dgraph);

    add_2q_num = 0;
    add_cnot_num = 0;
    add_swap_num = 0;
    add_bridge_num = 0;
    fidelity = 0;
    node_id = dgraph.nodeset.size();

    if (param_alpha == 0.0) param_alpha = 0.5;

    list<int> front_list;
    list<int> act_list;
    list<int> singlequbit_list;
    vector<bool> frozen(nqubits, false);

    FinalCircuit.nodeset.clear();

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    bool in_process = true;
    unsigned long long safety_counter = 0;

    // 다음 게이트가 inter인지 확인
    auto get_next_is_inter = [&](int logical_qubit) -> bool {
        auto& dl = Dlist[logical_qubit];
        if (dl.size() < 2) return false;
        auto it = next(dl.begin());
        int next_gate = *it;
        int nc = dgraph.nodeset[next_gate].control;
        int nt = dgraph.nodeset[next_gate].target;
        if (nc < 0 || nt < 0) return false;
        return extract_qpu_idx(nc) != extract_qpu_idx(nt);
    };

    // -------------------------------------------------------
    // MAIN LOOP
    // -------------------------------------------------------
    do {
        // === STEP 1: 즉시 실행 가능한 게이트 처리 ===
        bool progress = true;
        int inner_safety = 0;
        do {
            update_front_n_act_list(front_list, act_list, frozen);
            progress = check_direct_act_list(act_list, singlequbit_list,
                                             frozen, dgraph);
            act_list.unique();
            if (++inner_safety > 10000) break;
        } while (progress);

        // === STEP 2: act_list 전체 평가하여 최적 SWAP 선택 ===
        if (!act_list.empty())
        {
            struct Candidate {
                pair<int,int> swap;
                int gateid;
                double cost;
                bool is_inter;
            };
            vector<Candidate> all_candi;

            for (auto& gateid : act_list)
            {
                int control = dgraph.nodeset[gateid].control;
                int target  = dgraph.nodeset[gateid].target;
                int idx1 = extract_qpu_idx(control);
                int idx2 = extract_qpu_idx(target);

                if (layout_L[idx1].count(control) == 0) continue;
                if (layout_L[idx2].count(target)  == 0) continue;

                int Qc = layout_L[idx1][control];
                int Qt = layout_L[idx2][target];
                bool is_inter = (idx1 != idx2);

                bool next_is_inter_ctrl = get_next_is_inter(control);
                bool next_is_inter_tgt  = get_next_is_inter(target);

                vector<pair<int,int>> candi;

                if (is_inter) {
                    bool c_top = is_on_top_row(Qc, num_qubits, n);
                    bool t_top = is_on_top_row(Qt, num_qubits, n);
                    if (!c_top)
                        generate_inter_candidates(Qc, num_qubits, n,
                                                  multi_qpu_graph, qubit_Q, candi);
                    if (!t_top)
                        generate_inter_candidates(Qt, num_qubits, n,
                                                  multi_qpu_graph, qubit_Q, candi);
                }
                else {
                    generate_intra_candidates(Qc, Qt, idx1, idx2,
                                             control, target, candi, false);
                    if (candi.empty())
                        generate_intra_candidates(Qc, Qt, idx1, idx2,
                                                  control, target, candi, true);
                }

                for (auto& sw : candi)
                {
                    double cost = mapping_machine(is_inter, sw, dgraph, gateid);

                    // [개선 1] top row 큐비트를 끌어내리는 intra SWAP 패널티
                    if (!is_inter) {
                        auto top_eviction_penalty = [&](int Q_src, int Q_dst) {
                            if (is_on_top_row(Q_src, num_qubits, n) &&
                                !is_on_top_row(Q_dst, num_qubits, n))
                                cost -= 2.0;
                        };
                        top_eviction_penalty(sw.first, sw.second);
                        top_eviction_penalty(sw.second, sw.first);
                    }

                    // [개선 2] 다음 게이트가 inter면 top row 방향 보너스
                    if (!is_inter && (next_is_inter_ctrl || next_is_inter_tgt))
                    {
                        auto top_bonus = [&](int Q_focus) {
                            int Q_other = (sw.first == Q_focus) ? sw.second : sw.first;
                            if (!is_on_top_row(Q_focus, num_qubits, n) &&
                                is_on_top_row(Q_other, num_qubits, n))
                                cost += 1.5;
                            else {
                                int dist_before = dist_to_top(Q_focus, num_qubits, n);
                                int dist_after  = dist_to_top(Q_other, num_qubits, n);
                                if (dist_after < dist_before)
                                    cost += 1.0;
                            }
                        };
                        if (next_is_inter_ctrl) top_bonus(Qc);
                        if (next_is_inter_tgt)  top_bonus(Qt);
                    }

                    all_candi.push_back({sw, gateid, cost, is_inter});
                }
            }

            if (!all_candi.empty())
            {
                auto best_it = max_element(
                    all_candi.begin(), all_candi.end(),
                    [](const Candidate& a, const Candidate& b) {
                        return a.cost < b.cost;
                    });

                // 양수든 음수든 무조건 best 적용
                // (후보가 항상 존재하므로 진전 보장)
                do_swap(best_it->swap.first, best_it->swap.second);

                // cost 음수면 frozen 해제 + 막힌 게이트 뒤로
                if (best_it->cost <= 0) {
                    if (!act_list.empty()) {
                        int stuck_gate = act_list.front();
                        act_list.pop_front();
                        act_list.push_back(stuck_gate);
                    }
                    fill(frozen.begin(), frozen.end(), false);
                }
            }
            else {
                // 후보 자체가 없으면 막힌 게이트 뒤로 + frozen 해제
                if (!act_list.empty()) {
                    int stuck_gate = act_list.front();
                    act_list.pop_front();
                    act_list.push_back(stuck_gate);
                }
                fill(frozen.begin(), frozen.end(), false);
            }
        }
        else {
            fill(frozen.begin(), frozen.end(), false);
        }

        // === STEP 3: 종료 조건 ===
        int done_qubits = 0;
        for (int q = 0; q < (int)nqubits; q++)
            if (Dlist[q].empty()) done_qubits++;
        if (done_qubits == (int)nqubits)
            in_process = false;

        // === 안전장치 (진짜 버그 대비용으로만 유지) ===
        if (++safety_counter > 500000) {
            cout << "[WARN] safety counter triggered, breaking loop\n";
            add_2q_num = INT_MAX;
            break;
        }

    } while (in_process);

    cout << "[END] main_mapping completed. Added SWAP count: "
         << add_2q_num << "\n";
}

// =====================================================================
// 보조 함수들 (친구 코드에서 유지, 일부 보강)
// =====================================================================
void Qcircuit::QMapper::find_singlequbit_list(
    int gateid, list<int>& singlequbit_list, Circuit& dgraph)
{
    singlequbit_list.clear();
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;

    int erase_control = 0;
    int erase_target = 0;

    for (auto& id : Dlist_all[control]) {
        if (id == gateid) { erase_control++; break; }
        singlequbit_list.push_back(id);
        erase_control++;
    }
    for (int i = 0; i < erase_control; i++) Dlist_all[control].pop_front();

    for (auto& id : Dlist_all[target]) {
        if (id == gateid) { erase_target++; break; }
        singlequbit_list.push_back(id);
        erase_target++;
    }
    for (int i = 0; i < erase_target; i++) Dlist_all[target].pop_front();

    singlequbit_list.sort();
}

void Qcircuit::QMapper::update_act_dist2_list(
    list<int>& act_dist2_list, list<int>& act_list, Circuit& dgraph)
{
    act_dist2_list.clear();
    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);
        if (idx1 != idx2) continue;

        if (layout_L[idx1].count(control) == 0) continue;
        if (layout_L[idx2].count(target)  == 0) continue;

        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];
        if (multi_qpu_graph.dist[Q_control][Q_target] == 2)
            act_dist2_list.push_back(gateid);
    }
}

void Qcircuit::QMapper::move_qubit_up(int Q)
{
    // 호환성을 위한 더미 — 새 구조에서는 inter_candidates가 대체
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    int qpu_idx = Q / num_qubits;
    int local = Q % num_qubits;
    int row = local / n;
    int col = local % n;
    if (row == 0) return;

    int upper_local = (row - 1) * n + col;
    int upper_global = qpu_idx * num_qubits + upper_local;
    do_swap(Q, upper_global);
}

void Qcircuit::QMapper::generate_candi_list(
    list<int>& act_list,
    vector<pair<pair<int,int>, int>>& candi_list,
    Circuit& dgraph)
{
    // 호환성을 위한 더미 — 새 구조에서는 main_mapping 내부에서 직접 처리
    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);
        if (idx1 != idx2) continue;  // intra only

        if (layout_L[idx1].count(control) == 0) continue;
        if (layout_L[idx2].count(target)  == 0) continue;

        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];

        for (int i = 0; i < multi_qpu_graph.node_size; i++)
        {
            if (extract_qpu_idx(i) != idx1) continue;
            if (multi_qpu_graph.dist[Q_control][i] == 1) {
                if (cal_SWAP_effect(control, target, i, Q_control) > 0)
                    candi_list.push_back({{min(i, Q_control), max(i, Q_control)}, gateid});
            }
            if (multi_qpu_graph.dist[Q_target][i] == 1) {
                if (cal_SWAP_effect(control, target, i, Q_target) > 0)
                    candi_list.push_back({{min(i, Q_target), max(i, Q_target)}, gateid});
            }
        }
    }
}

bool Qcircuit::QMapper::top_direction_swap_filter(int Q_from, int Q_to)
{
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    return dist_to_top(Q_to, num_qubits, n) <
           dist_to_top(Q_from, num_qubits, n);
}