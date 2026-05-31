#include "circuit.h"
#include <algorithm>
#include <iostream>
#include <cmath>
#include <climits>

using namespace std;
using namespace Qcircuit;

#define EMPTY_NODE -1

static const double P_2Q_ERROR   = 1e-3;                       // 2-qubit gate error
static const double F_LOCAL_2Q   = 1.0 - P_2Q_ERROR;           // 로컬 2q fidelity
static const double F_SWAP       = F_LOCAL_2Q * F_LOCAL_2Q * F_LOCAL_2Q; // SWAP=CNOT*3
static const double F_EPR        = 0.85;                       // inter-QPU 원격 게이트

// 비용 가중치
static const double WEIGHT_SWAP            = 3.0;
static const double WEIGHT_VIRTUAL_BUFFER  = 10.0;

// 헬퍼들
static inline int safe_get_logical(
    vector<map<int,int>>& qubit_Q, int qpu_idx, int physical)
{
    if (qpu_idx < 0 || qpu_idx >= (int)qubit_Q.size()) return -1;
    auto it = qubit_Q[qpu_idx].find(physical);
    if (it == qubit_Q[qpu_idx].end()) return -1;
    return it->second;
}

static inline bool is_on_top_row(int Q, int num_qubits, int n)
{
    int local = Q % num_qubits;
    return (local / n) == 0;
}

static inline int dist_to_top(int Q, int num_qubits, int n)
{
    int local = Q % num_qubits;
    return local / n;
}

// mapping_machine - 통합 비용 함수
double Qcircuit::QMapper::mapping_machine(
    bool is_inter_gate, pair<int,int> SWAP_pair,
    Circuit& dgraph, int gateid)
{
    int Q1 = SWAP_pair.first;
    int Q2 = SWAP_pair.second;
    int qpu_idx1 = Q1 / num_qubits;
    int qpu_idx2 = Q2 / num_qubits;

    int q1 = safe_get_logical(qubit_Q, qpu_idx1, Q1);
    int q2 = safe_get_logical(qubit_Q, qpu_idx2, Q2);

    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    if (is_inter_gate)
    {
        double final_cost = 0.0;
        if (q1 != -1 && (q1 == control || q1 == target)) {
            int dist_before = dist_to_top(Q1, num_qubits, n);
            int dist_after  = dist_to_top(Q2, num_qubits, n);
            double current_cost = (dist_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost  = (dist_after  * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            final_cost += (current_cost - future_cost);
        }
        if (q2 != -1 && (q2 == control || q2 == target)) {
            int dist_before = dist_to_top(Q2, num_qubits, n);
            int dist_after  = dist_to_top(Q1, num_qubits, n);
            double current_cost = (dist_before * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            double future_cost  = (dist_after  * WEIGHT_SWAP) + WEIGHT_VIRTUAL_BUFFER;
            final_cost += (current_cost - future_cost);
        }
        if ((q1 == -1 && q2 != -1 && (q2 == control || q2 == target)) ||
            (q2 == -1 && q1 != -1 && (q1 == control || q1 == target))) {
            if (final_cost == 0.0) final_cost = 1.0;
        }
        // inter SWAP의 경우 dist_to_top 비용만 사용 (단순화)
        return final_cost;
    }
    else
    {
        if (q1 == -1 || q2 == -1) return 0.5;
        return cal_MCPE(SWAP_pair, dgraph);
    }
}

void Qcircuit::QMapper::update_front_n_act_list(
    list<int>& front_list, list<int>& act_list, vector<bool>& frozen)
{
    for (int q = 0; q < (int)nqubits; q++)
    {
        if (frozen[q]) continue;
        if (Dlist[q].empty()) continue;
        int gateid = Dlist[q].front();
        auto front_it = find(front_list.begin(), front_list.end(), gateid);
        if (front_it != front_list.end()) {
            front_list.erase(front_it);
            if (find(act_list.begin(), act_list.end(), gateid) == act_list.end())
                act_list.push_back(gateid);
        }
        else front_list.push_back(gateid);
        frozen[q] = true;
    }
}

//   버퍼 n개 제약 추가
bool Qcircuit::QMapper::check_direct_act_list(
    list<int>& act_list, list<int>& singlequbit_list,
    vector<bool>& frozen, Circuit& dgraph)
{
    bool complete = false;
    vector<int> erase_list;
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));

    // 이번 라운드에서 각 QPU의 버퍼(top row)를 몇 개 점유했는지
    vector<int> buffer_usage(num_qpu, 0);

    for (auto& gateid : act_list)
    {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);

        if (idx1 < 0 || idx2 < 0) continue;  // 매핑 안 된 큐비트 스킵
        if (idx1 >= (int)layout_L.size() || idx2 >= (int)layout_L.size()) continue;
        if (layout_L[idx1].count(control) == 0) continue;
        if (layout_L[idx2].count(target)  == 0) continue;

        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];

        if (Q_control < 0 || Q_control >= multi_qpu_graph.node_size) continue;
        if (Q_target  < 0 || Q_target  >= multi_qpu_graph.node_size) continue;

        bool can_execute = false;

        if (idx1 == idx2) {
            // Intra: 거리 1이면 실행 (버퍼 제약 무관)
            if (multi_qpu_graph.dist[Q_control][Q_target] == 1)
                can_execute = true;
        }
        else {
            // Inter: 양쪽 모두 top row + 버퍼 여유 있어야 실행
            bool c_top = is_on_top_row(Q_control, num_qubits, n);
            bool t_top = is_on_top_row(Q_target,  num_qubits, n);

            if (c_top && t_top) {
                // 버퍼 용량 검사 (각 QPU는 n개 버퍼슬롯)
                if (buffer_usage[idx1] < n && buffer_usage[idx2] < n) {
                    can_execute = true;
                    buffer_usage[idx1]++;
                    buffer_usage[idx2]++;
                }
                // 용량 초과 시 이번 라운드는 실행 보류
            }
        }

        if (can_execute) {
            erase_list.push_back(gateid);
            FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
            if (!Dlist[control].empty()) Dlist[control].pop_front();
            if (!Dlist[target].empty())  Dlist[target].pop_front();
            frozen[control] = false;
            frozen[target]  = false;
            complete = true;

            //  inter-gate 실행 카운트 (add_bridge_num 재활용)
            if (idx1 != idx2) add_bridge_num++;
        }
    }

    for (auto gid : erase_list)
        act_list.remove(gid);

    return complete;
}

void Qcircuit::QMapper::do_swap(int Q1, int Q2)
{
    // [안정성] 물리 큐비트 범위 가드
    if (Q1 < 0 || Q1 >= multi_qpu_graph.node_size) return;
    if (Q2 < 0 || Q2 >= multi_qpu_graph.node_size) return;

    int qpu_idx1 = Q1 / num_qubits;
    int qpu_idx2 = Q2 / num_qubits;
    if (qpu_idx1 < 0 || qpu_idx1 >= (int)qubit_Q.size()) return;
    if (qpu_idx2 < 0 || qpu_idx2 >= (int)qubit_Q.size()) return;

    int q1 = safe_get_logical(qubit_Q, qpu_idx1, Q1);
    int q2 = safe_get_logical(qubit_Q, qpu_idx2, Q2);

    if (q1 != -1) {
        int core1 = extract_qpu_idx(q1);
        if (core1 >= 0 && core1 < (int)layout_L.size()) layout_L[core1][q1] = Q2;
    }
    if (q2 != -1) {
        int core2 = extract_qpu_idx(q2);
        if (core2 >= 0 && core2 < (int)layout_L.size()) layout_L[core2][q2] = Q1;
    }

    if (q2 != -1) qubit_Q[qpu_idx1][Q1] = q2;
    else          qubit_Q[qpu_idx1].erase(Q1);
    if (q1 != -1) qubit_Q[qpu_idx2][Q2] = q1;
    else          qubit_Q[qpu_idx2].erase(Q2);

    if (q1 != -1 || q2 != -1) {
        add_swap(Q1, Q2, FinalCircuit);
        add_2q_num++;
    }
}

void Qcircuit::QMapper::generate_intra_candidates(
    int Q_control, int Q_target, int idx1, int idx2,
    int control, int target,
    vector<pair<int,int>>& candi, bool force_mode)
{
    for (int i = 0; i < multi_qpu_graph.node_size; i++)
    {
        if (multi_qpu_graph.dist[Q_control][i] == 1)
        {
            if (!force_mode && extract_qpu_idx(i) != idx1 && (i / num_qubits) != idx1) {
            } else {
                int i_qpu = i / num_qubits;
                bool is_empty = (i_qpu >= 0 && i_qpu < (int)qubit_Q.size() &&
                                 qubit_Q[i_qpu].count(i) == 0);
                if (is_empty)
                    candi.push_back({min(i, Q_control), max(i, Q_control)});
                else if (force_mode ||
                         cal_SWAP_effect(control, target, i, Q_control) > 0)
                    candi.push_back({min(i, Q_control), max(i, Q_control)});
            }
        }
        if (multi_qpu_graph.dist[Q_target][i] == 1)
        {
            if (!force_mode && extract_qpu_idx(i) != idx2 && (i / num_qubits) != idx2) {
            } else {
                int i_qpu = i / num_qubits;
                bool is_empty = (i_qpu >= 0 && i_qpu < (int)qubit_Q.size() &&
                                 qubit_Q[i_qpu].count(i) == 0);
                if (is_empty)
                    candi.push_back({min(i, Q_target), max(i, Q_target)});
                else if (force_mode ||
                         cal_SWAP_effect(control, target, i, Q_target) > 0)
                    candi.push_back({min(i, Q_target), max(i, Q_target)});
            }
        }
    }
    sort(candi.begin(), candi.end());
    candi.erase(unique(candi.begin(), candi.end()), candi.end());
}

static void generate_inter_candidates(
    int Q_focus, int num_qubits, int n,
    Graph& multi_qpu_graph, vector<map<int,int>>& qubit_Q,
    vector<pair<int,int>>& candi)
{
    int qpu_idx = Q_focus / num_qubits;
    int local = Q_focus % num_qubits;
    int row = local / n;
    int col = local % n;
    if (row == 0) return;

    vector<int> neighbors;
    int up = qpu_idx * num_qubits + (row - 1) * n + col;
    neighbors.push_back(up);
    if (col > 0) {
        int up_left = qpu_idx * num_qubits + (row - 1) * n + (col - 1);
        if (multi_qpu_graph.dist[Q_focus][up_left] == 1)
            neighbors.push_back(up_left);
    }
    if (col < n - 1) {
        int up_right = qpu_idx * num_qubits + (row - 1) * n + (col + 1);
        if (multi_qpu_graph.dist[Q_focus][up_right] == 1)
            neighbors.push_back(up_right);
    }
    for (int pos : neighbors) {
        if (pos / num_qubits != qpu_idx) continue;
        if (multi_qpu_graph.dist[Q_focus][pos] != 1) continue;
        candi.push_back({min(Q_focus, pos), max(Q_focus, pos)});
    }
}

bool Qcircuit::QMapper::handle_intra_gate(
    int gateid, Circuit& dgraph,
    vector<bool>& frozen, list<int>& act_list)
{
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;
    int idx1 = extract_qpu_idx(control);
    int idx2 = extract_qpu_idx(target);
    if (idx1 < 0 || idx2 < 0) return false;
    if (idx1 >= (int)layout_L.size() || idx2 >= (int)layout_L.size()) return false;
    if (layout_L[idx1].count(control) == 0) return false;
    if (layout_L[idx2].count(target)  == 0) return false;
    int Q_control = layout_L[idx1][control];
    int Q_target  = layout_L[idx2][target];

    if (multi_qpu_graph.dist[Q_control][Q_target] == 1) {
        FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
        if (control < (int)Dlist.size() && !Dlist[control].empty()) Dlist[control].pop_front();
        if (target  < (int)Dlist.size() && !Dlist[target].empty())  Dlist[target].pop_front();
        frozen[control] = false;
        frozen[target]  = false;
        act_list.remove(gateid);
        return true;
    }
    vector<pair<int,int>> candi;
    generate_intra_candidates(Q_control, Q_target, idx1, idx2,
                              control, target, candi, false);
    if (candi.empty())
        generate_intra_candidates(Q_control, Q_target, idx1, idx2,
                                  control, target, candi, true);
    if (candi.empty()) return false;
    double best_cost = -1e18;
    pair<int,int> best_swap = candi[0];
    for (auto& sw : candi) {
        double cost = mapping_machine(false, sw, dgraph, gateid);
        if (cost > best_cost) { best_cost = cost; best_swap = sw; }
    }
    do_swap(best_swap.first, best_swap.second);
    return false;
}

bool Qcircuit::QMapper::handle_inter_gate(
    int gateid, Circuit& dgraph,
    vector<bool>& frozen, list<int>& act_list)
{
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;
    int idx1 = extract_qpu_idx(control);
    int idx2 = extract_qpu_idx(target);
    if (idx1 < 0 || idx2 < 0) return false;
    if (idx1 >= (int)layout_L.size() || idx2 >= (int)layout_L.size()) return false;
    if (layout_L[idx1].count(control) == 0) return false;
    if (layout_L[idx2].count(target)  == 0) return false;
    int Qc = layout_L[idx1][control];
    int Qt = layout_L[idx2][target];
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    bool c_top = is_on_top_row(Qc, num_qubits, n);
    bool t_top = is_on_top_row(Qt, num_qubits, n);

    if (c_top && t_top) {
        FinalCircuit.nodeset.push_back(dgraph.nodeset[gateid]);
        if (control < (int)Dlist.size() && !Dlist[control].empty()) Dlist[control].pop_front();
        if (target  < (int)Dlist.size() && !Dlist[target].empty())  Dlist[target].pop_front();
        frozen[control] = false;
        frozen[target]  = false;
        act_list.remove(gateid);
        add_bridge_num++;  // [ADD-2] inter-gate 실행 카운트
        return true;
    }
    vector<int> need_up;
    if (!c_top) need_up.push_back(Qc);
    if (!t_top) need_up.push_back(Qt);
    vector<pair<int,int>> candi;
    for (int Q : need_up)
        generate_inter_candidates(Q, num_qubits, n,
                                  multi_qpu_graph, qubit_Q, candi);
    sort(candi.begin(), candi.end());
    candi.erase(unique(candi.begin(), candi.end()), candi.end());
    if (candi.empty()) return false;
    double best_cost = -1e18;
    pair<int,int> best_swap = candi[0];
    for (auto& sw : candi) {
        double cost = mapping_machine(true, sw, dgraph, gateid);
        if (cost > best_cost) { best_cost = cost; best_swap = sw; }
    }
    if (best_cost <= 0) return false;
    do_swap(best_swap.first, best_swap.second);
    return false;
}

static int compute_fidelity_percent(
    int add_swap_count,
    int inter_exec_count)
{
    // 안전: 음수 방지
    if (add_swap_count  < 0) add_swap_count  = 0;
    if (inter_exec_count < 0) inter_exec_count = 0;

    // log 공간 계산 (언더플로 방지)
    double log_fid = 0.0;
    log_fid += add_swap_count   * std::log(F_SWAP);
    log_fid += inter_exec_count * std::log(F_EPR);

    double F = std::exp(log_fid);

    // 0~1 clamp
    if (F < 0.0) F = 0.0;
    if (F > 1.0) F = 1.0;

    // 퍼센트 정수로 반올림 (run.cpp extract_int가 소수점 못 읽으므로)
    int percent = static_cast<int>(std::round(F * 100.0));

    cout << "[FIDELITY] swaps=" << add_swap_count
         << " inter_epr=" << inter_exec_count
         << " => F=" << F << " (" << percent << "%)\n";

    return percent;
}

void Qcircuit::QMapper::main_mapping(
    Circuit& dgraph, bool BRIDGE_MODE, int n_buffer)
{
    cout << "[START] main_mapping v5\n";

    make_Dlist(dgraph);
    make_Dlist_all(dgraph);

    add_2q_num = 0;
    add_cnot_num = 0;
    add_swap_num = 0;
    add_bridge_num = 0;
    fidelity = 0;
    node_id = dgraph.nodeset.size();
    if (param_alpha == 0.0) param_alpha = 0.6;  // [개선] 0.5→0.7: 미래 게이트 영향 더 반영

    list<int> front_list;
    list<int> act_list;
    list<int> singlequbit_list;
    vector<bool> frozen(nqubits, false);
    FinalCircuit.nodeset.clear();

    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    bool in_process = true;
    unsigned long long safety_counter = 0;

    // [개선] 다음 k개 게이트 중 inter가 있으면 true (기존: 다음 1개만 확인)
    auto get_next_is_inter = [&](int logical_qubit) -> bool {
        auto& dl = Dlist[logical_qubit];
        if (dl.size() < 2) return false;
        int lookahead = 3;  // 다음 3개 게이트까지 확인
        auto it = next(dl.begin());
        for (int k = 0; k < lookahead && it != dl.end(); ++k, ++it) {
            int ng = *it;
            int nc = dgraph.nodeset[ng].control;
            int nt = dgraph.nodeset[ng].target;
            if (nc < 0 || nt < 0) continue;
            if (extract_qpu_idx(nc) != extract_qpu_idx(nt)) return true;
        }
        return false;
    };

    do {
        // STEP 1
        bool progress = true;
        int inner_safety = 0;
        do {
            update_front_n_act_list(front_list, act_list, frozen);
            progress = check_direct_act_list(act_list, singlequbit_list,
                                             frozen, dgraph);
            act_list.unique();
            if (++inner_safety > 10000) break;
        } while (progress);

        // STEP 2
        if (!act_list.empty())
        {
            struct Candidate {
                pair<int,int> swap;
                int gateid;
                double cost;
                bool is_inter;
            };
            vector<Candidate> all_candi;

            // 이번 라운드 버퍼 점유 예측 (inter SWAP 후보 평가에 반영)
            vector<int> buffer_pressure(num_qpu, 0);
            for (auto& gateid : act_list) {
                int c = dgraph.nodeset[gateid].control;
                int t = dgraph.nodeset[gateid].target;
                int ic = extract_qpu_idx(c);
                int it = extract_qpu_idx(t);
                if (ic >= 0 && it >= 0 && ic != it) {
                    buffer_pressure[ic]++;
                    buffer_pressure[it]++;
                }
            }

            for (auto& gateid : act_list)
            {
                int control = dgraph.nodeset[gateid].control;
                int target  = dgraph.nodeset[gateid].target;
                int idx1 = extract_qpu_idx(control);
                int idx2 = extract_qpu_idx(target);
                if (idx1 < 0 || idx2 < 0) continue;
                if (idx1 >= (int)layout_L.size() || idx2 >= (int)layout_L.size()) continue;
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

                    // [개선] top eviction penalty: inter-gate를 가진 큐비트는 강하게 보호
                    // inter-gate 없는 큐비트는 부드러운 패널티 (radd_250 같이 inter 많은 회로 대응)
                    if (!is_inter) {
                        auto top_eviction_penalty = [&](int Q_src, int Q_dst) {
                            if (is_on_top_row(Q_src, num_qubits, n) &&
                                !is_on_top_row(Q_dst, num_qubits, n)) {
                                int qpu_s = Q_src / num_qubits;
                                int lq = safe_get_logical(qubit_Q, qpu_s, Q_src);
                                bool has_inter = false;
                                if (lq >= 0 && lq < (int)Dlist.size() && !Dlist[lq].empty()) {
                                    int ng = Dlist[lq].front();
                                    int nc = dgraph.nodeset[ng].control;
                                    int nt = dgraph.nodeset[ng].target;
                                    if (nc >= 0 && nt >= 0 && extract_qpu_idx(nc) != extract_qpu_idx(nt))
                                        has_inter = true;
                                }
                                // inter gate 가진 큐비트: 강한 보호 / 없으면 약한 패널티
                                cost -= has_inter ? 4.0 : 1.5;
                            }
                        };
                        top_eviction_penalty(sw.first, sw.second);
                        top_eviction_penalty(sw.second, sw.first);
                    }

                    // [개선] next-gate inter bonus 강화: 1.5→3.0, 1.0→2.0
                    if (!is_inter && (next_is_inter_ctrl || next_is_inter_tgt)) {
                        auto top_bonus = [&](int Q_focus) {
                            int Q_other = (sw.first == Q_focus) ? sw.second : sw.first;
                            if (!is_on_top_row(Q_focus, num_qubits, n) &&
                                is_on_top_row(Q_other, num_qubits, n))
                                cost += 2.0;
                            else {
                                int db = dist_to_top(Q_focus, num_qubits, n);
                                int da = dist_to_top(Q_other, num_qubits, n);
                                if (da < db) cost += 1.0;
                            }
                        };
                        if (next_is_inter_ctrl) top_bonus(Qc);
                        if (next_is_inter_tgt)  top_bonus(Qt);
                    }

                    if (is_inter) {
                        if (idx1 >= 0 && buffer_pressure[idx1] > n)
                            cost -= 0.5 * (buffer_pressure[idx1] - n);
                        if (idx2 >= 0 && buffer_pressure[idx2] > n)
                            cost -= 0.5 * (buffer_pressure[idx2] - n);
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

                do_swap(best_it->swap.first, best_it->swap.second);

                if (best_it->cost <= 0) {
                    // [개선] 전체 frozen clear 대신 stuck gate의 큐비트만 해제
                    if (!act_list.empty()) {
                        int stuck = act_list.front();
                        act_list.pop_front();
                        act_list.push_back(stuck);
                        // stuck gate 관련 큐비트만 unfreeze
                        int sc = dgraph.nodeset[stuck].control;
                        int st = dgraph.nodeset[stuck].target;
                        if (sc >= 0 && sc < (int)frozen.size()) frozen[sc] = false;
                        if (st >= 0 && st < (int)frozen.size()) frozen[st] = false;
                    } else {
                        fill(frozen.begin(), frozen.end(), false);
                    }
                }
            }
            else {
                if (!act_list.empty()) {
                    int stuck = act_list.front();
                    act_list.pop_front();
                    act_list.push_back(stuck);
                }
                fill(frozen.begin(), frozen.end(), false);
            }
        }
        else {
            fill(frozen.begin(), frozen.end(), false);
        }

        // 종료 조건
        int done_qubits = 0;
        for (int q = 0; q < (int)nqubits; q++)
            if (Dlist[q].empty()) done_qubits++;
        if (done_qubits == (int)nqubits)
            in_process = false;

        if (++safety_counter > 2000000ULL) {
            cout << "[WARN] safety counter triggered\n";
            add_2q_num = INT_MAX;
            break;
        }
    } while (in_process);

    //   Fidelity 계산
    //   add_bridge_num = 실행된 inter-gate 수
    //   add_swap_num   = 추가된 SWAP 수
    if (add_2q_num != INT_MAX) {
        int fid_percent = compute_fidelity_percent(add_swap_num, add_bridge_num);
        fidelity = static_cast<float>(fid_percent);  // 정수 퍼센트 저장
    } else {
        fidelity = 0;  // OOB or 데드락은 0
    }

    cout << "[END] main_mapping v5.1 SWAPs=" << add_2q_num
         << " inter_exec=" << add_bridge_num
         << " fidelity=" << fidelity << "%\n";
}

void Qcircuit::QMapper::find_singlequbit_list(
    int gateid, list<int>& singlequbit_list, Circuit& dgraph)
{
    singlequbit_list.clear();
    int control = dgraph.nodeset[gateid].control;
    int target  = dgraph.nodeset[gateid].target;
    int erase_control = 0, erase_target = 0;
    for (auto& id : Dlist_all[control]) {
        if (id == gateid) { erase_control++; break; }
        singlequbit_list.push_back(id); erase_control++;
    }
    for (int i = 0; i < erase_control; i++) Dlist_all[control].pop_front();
    for (auto& id : Dlist_all[target]) {
        if (id == gateid) { erase_target++; break; }
        singlequbit_list.push_back(id); erase_target++;
    }
    for (int i = 0; i < erase_target; i++) Dlist_all[target].pop_front();
    singlequbit_list.sort();
}

void Qcircuit::QMapper::update_act_dist2_list(
    list<int>& act_dist2_list, list<int>& act_list, Circuit& dgraph)
{
    act_dist2_list.clear();
    for (auto& gateid : act_list) {
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
    for (auto& gateid : act_list) {
        int control = dgraph.nodeset[gateid].control;
        int target  = dgraph.nodeset[gateid].target;
        int idx1 = extract_qpu_idx(control);
        int idx2 = extract_qpu_idx(target);
        if (idx1 != idx2) continue;
        if (layout_L[idx1].count(control) == 0) continue;
        if (layout_L[idx2].count(target)  == 0) continue;
        int Q_control = layout_L[idx1][control];
        int Q_target  = layout_L[idx2][target];
        for (int i = 0; i < multi_qpu_graph.node_size; i++) {
            if (extract_qpu_idx(i) != idx1) continue;
            if (multi_qpu_graph.dist[Q_control][i] == 1 &&
                cal_SWAP_effect(control, target, i, Q_control) > 0)
                candi_list.push_back({{min(i, Q_control), max(i, Q_control)}, gateid});
            if (multi_qpu_graph.dist[Q_target][i] == 1 &&
                cal_SWAP_effect(control, target, i, Q_target) > 0)
                candi_list.push_back({{min(i, Q_target), max(i, Q_target)}, gateid});
        }
    }
}

bool Qcircuit::QMapper::top_direction_swap_filter(int Q_from, int Q_to)
{
    int n = static_cast<int>(std::floor(std::sqrt(num_qubits)));
    return dist_to_top(Q_to, num_qubits, n) <
           dist_to_top(Q_from, num_qubits, n);
}