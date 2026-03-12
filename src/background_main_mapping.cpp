//mappingFunction.cpp
#include "circuit.h"

using namespace std;
using namespace Qcircuit;

bool compare_by_cost(pair< pair<int, int>, pair<int, int> > p1, pair< pair<int, int>, pair<int, int> > p2)
{
    return p1.second.second > p2.second.second;
}

int Qcircuit::QMapper::cal_SWAP_effect(const int q1, const int q2, const int Q1, const int Q2)
{
    int idx1,idx2,idx3,idx4;
    idx1=extract_qpu_idx(q1);
    idx2=extract_qpu_idx(q2);
    idx3=Q1/num_qubits;
    idx4=Q2/num_qubits;
    int controlQ = layout_L[idx1][q1];
    int targetQ = layout_L[idx2][q2];

    int dist_before_swap = multi_qpu_graph.dist[controlQ][targetQ];
    int swapq1 = qubit_Q[idx3][Q1];
    int swapq2 = qubit_Q[idx4][Q2];

    if(controlQ == Q1 || controlQ == Q2)
        controlQ = (controlQ == Q1) ? Q2 : Q1;
    if(targetQ == Q1 || targetQ == Q2)
        targetQ = (targetQ == Q1) ? Q2 : Q1;

    int dist_after_swap  = multi_qpu_graph.dist[controlQ][targetQ];
    int effect = dist_before_swap - dist_after_swap;

    return effect;
}

double Qcircuit::QMapper::cal_MCPE(const pair<int, int> p, Circuit& dgraph)
{
    int idx1 = p.first/num_qubits;
    int idx2 = p.second/num_qubits;
    int q1 = qubit_Q[idx1][p.first];
    int q2 = qubit_Q[idx2][p.second];

    double MCPE = 0;
    int dist = 0;
    double power = 1.0;

    if(q1 < nqubits)
        for(auto& gateid : Dlist[q1])
        {
            Gate &g = dgraph.nodeset[gateid];
            int control = g.control;
            int target  = g.target;
            int effect  = cal_SWAP_effect(control, target, p.first, p.second);
            if(effect < 0) break;
            if(power < 0.0000001) break;
            MCPE += effect * power;
            power *= param_alpha;
            dist++;
        }
    
    dist = 0;
    power = 1.0;
    if(q2 < nqubits)
        for(auto& gateid : Dlist[q2])
        {
            Gate &g = dgraph.nodeset[gateid];
            int control = g.control;
            int target  = g.target;
            int effect  = cal_SWAP_effect(control, target, p.first, p.second);
            if(effect < 0) break;
            if(power < 0.0000001) break;
            MCPE += effect * power;
            power *= param_alpha;
            dist++;
        }

    return MCPE*10;
}

void Qcircuit::QMapper::find_max_cost(pair<int, int>& SWAP, int& gateid, double& max_cost,
                                             vector< pair< pair<int, int>, pair<int, double> > >& v, list<int>& act_list,
                                             vector< pair< pair<int, int>, pair<int, double> > >& history)
{
    sort(v.begin(), v.end(), compare_by_cost);
    
    int index = 0;
    pair< pair<int, int>, pair<int, double> >& ref = v[index];

    for(auto pair : history)
        if(pair == ref)
            index++;

    if(index >= v.size()) index = v.size()-1;

    SWAP = v[index].first;
    gateid = v[index].second.first;
    max_cost = v[index].second.second;
}

void Qcircuit::QMapper::layout_swap(const int b1, const int b2)
{
    int idx1 = extract_qpu_idx(b1);
    int idx2 = extract_qpu_idx(b2);
    int& q1 = layout_L[idx1][b1];
    int& q2 = layout_L[idx2][b2];

    
    if(idx1 == idx2)
    {
        swap(q1, q2);
        swap(qubit_Q[idx1][q1], qubit_Q[idx1][q2]);
    }
        
}

void Qcircuit::QMapper::add_cnot(int c_qubit, int t_qubit, Circuit& graph)
{
    Qcircuit::Gate cnot;
    
    cnot.id = node_id;
    cnot.control = c_qubit;
    cnot.target = t_qubit;
    cnot.type = GATETYPE::CNOT;
    
    graph.nodeset.push_back(cnot);
    add_cnot_num++;
    node_id++;
}

void Qcircuit::QMapper::add_swap(int b1, int b2, Circuit& graph)
{
    add_cnot(b1, b2, graph);
    add_cnot(b2, b1, graph);
    add_cnot(b1, b2, graph);

    add_swap_num++;
}

void Qcircuit::QMapper::add_bridge(int qs, int qt, Circuit& graph)
{
    int qb;
    int bs, bt, b;
    /*qs: control, qt: target*/

    for(int i=0; i<positions; i++)
    {
        if(multi_qpu_graph.dist[qs][i] == 1)
        {
            if(multi_qpu_graph.dist[qt][i] == 1)
                qb = i;
            else continue;
        }
        else continue;
    }
    int idx1 = qs/num_qubits;
    int idx2 = qt/num_qubits;
    int idx3 = qb/num_qubits;

    bs = qubit_Q[idx1][qs];
    bt = qubit_Q[idx2][qt];
    b = qubit_Q[idx3][qb];
    
    add_cnot(b, bt, graph);
    add_cnot(bs, b, graph);
    add_cnot(b, bt, graph);
    add_cnot(bs, b, graph);
    
    add_cnot_num--;
    add_bridge_num++;
}