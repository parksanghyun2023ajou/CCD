#include "circuit.h"
using namespace std;
using namespace Qcircuit;

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

void Qcircuit::QMapper::make_CNOT(bool i)
{
    Dgraph_cnot.nodeset.clear();

    if(i)
    {
        for(const auto& gate : Dgraph.nodeset)
        {
            if(gate.type != GATETYPE::CNOT) continue;
            const int target  = gate.target;
            const int control = gate.control;
            add_cnot(control, target, Dgraph_cnot);
            add_cnot_num--;
        }
    }
}

Graph Qcircuit::QMapper::make_interactionGraph(bool i)
{
    bool** gen_edge = new bool*[nqubits]; //edge generated (for interaction graph)
    for(int i = 0; i < nqubits; i++)
        gen_edge[i] = new bool[nqubits];

    //gen_edge initialize
    for(int i = 0; i < nqubits; i++)
        for(int j = 0; j < nqubits; j++)
            if(i == j)
                gen_edge[i][j] = true;
            else
                gen_edge[i][j] = false;
    Graph g;
    for(int i = 0; i < nqubits; i++)
        g.addnode();
    for(const auto& gate : Dgraph_cnot.nodeset)
    {
        if(gate.type != GATETYPE::CNOT) continue;
        const int id = gate.id;
        const int target = gate.target;
        const int control = gate.control;
        if(gen_edge[target][control]) continue;
        g.addedge(target, control, id);
        gen_edge[target][control] = true;
        gen_edge[control][target] = true;
    }
    
    for(int i = 0; i < nqubits; i++)
        delete[] gen_edge[i];
    delete[] gen_edge;

    return g;
}

Graph Qcircuit::QMapper::make_interactionNumberGraph(bool i)
{
    bool** gen_edge = new bool*[nqubits]; //edge generated (for interaction graph)
    for(int i = 0; i < nqubits; i++)
        gen_edge[i] = new bool[nqubits];

    //gen_edge initialize
    for(int i = 0; i < nqubits; i++)
        for(int j = 0; j < nqubits; j++)
            if(i == j)
                gen_edge[i][j] = true;
            else
                gen_edge[i][j] = false;

    Graph g;
    for(int i = 0; i < nqubits; i++)
        g.addnode();
    for(const auto& gate : Dgraph_cnot.nodeset)
    {
        if(gate.type != GATETYPE::CNOT) continue;
        const int id = gate.id;
        const int target = gate.target;
        const int control = gate.control;
        if(gen_edge[target][control])
        {
            for(auto& e : g.edgeset)
            {
                if(e.second.gettargetid()==target && e.second.getsourceid()==control)
                {
                    int weight = e.second.getweight();
                    e.second.setweight(weight+1);
                }
                else if(e.second.gettargetid()==control && e.second.getsourceid()==target)
                {
                    int weight = e.second.getweight();
                    e.second.setweight(weight+1);
                }

                else continue;
            }
        }
        else
        g.addedge(target, control, 1);
        gen_edge[target][control] = true;
        gen_edge[control][target] = true;

    }
    
    for(int i = 0; i < nqubits; i++)
        delete[] gen_edge[i];
    delete[] gen_edge;

    return g;
}


Graph Qcircuit::QMapper::make_interactionMixgraph(bool i, int n)
{
    bool** gen_edge = new bool*[nqubits]; //edge generated (for interaction graph)
    for(int i = 0; i < nqubits; i++)
        gen_edge[i] = new bool[nqubits];

    //gen_edge initialize
    for(int i = 0; i < nqubits; i++)
        for(int j = 0; j < nqubits; j++)
            if(i == j)
                gen_edge[i][j] = true;
            else
                gen_edge[i][j] = false;

    Graph g;
    for(int i = 0; i < nqubits; i++)
        g.addnode();

    int dgraphSize = Dgraph_cnot.nodeset.size();
    int m = dgraphSize / n;
    int remainder = dgraphSize % n;
    
    vector<double> costWeightTable;
    double costParam = 0.9;
    double cost = 1.0;
    costWeightTable.clear();
    for(int i=0; i<n+1; i++)
    {
        costWeightTable.push_back(cost);
        cost *= costParam;
    }
    
    
    int for_i=0;
    int for_n=0;
    double weight_cost;
    for(const auto& gate : Dgraph_cnot.nodeset)
    {
        if(for_i==m)
        {
            for_i=0;
            for_n++;
        }
        weight_cost = costWeightTable[for_n];
        const int id = gate.id;
        const int target = gate.target;
        const int control = gate.control;
        if(gen_edge[target][control])
        {
            for(auto& e : g.edgeset)
            {
                if(e.second.gettargetid()==target && e.second.getsourceid()==control)
                {
                    double weight = e.second.getweight();
                    e.second.setweight(weight+weight_cost*10);
                }
                else if(e.second.gettargetid()==control && e.second.getsourceid()==target)
                {
                    double weight = e.second.getweight();
                    e.second.setweight(weight+weight_cost*10);
                }

                else continue;
            }
        }
        else
        g.addedge(target, control, weight_cost);
        gen_edge[target][control] = true;
        gen_edge[control][target] = true;

        for_i++;
    }

    for(int i = 0; i < nqubits; i++)
        delete[] gen_edge[i];
    delete[] gen_edge;

    return g;

}

void Qcircuit::QMapper::print_interactionGraph( Graph& g)
{
    cout << "\n================ Interaction Graph ================\n";

    cout << "[Nodes]\n";
    for ( auto& [id, node] : g.nodeset)
    {
        cout << "  Node " << id << " (weight=" << node.getweight() << "): edges -> ";
        for (int eid : node.edges)
            cout << eid << " ";
        cout << "\n";
    }

    cout << "\n[Edges]\n";
    for (const auto& [eid, edge] : g.edgeset)
    {
        cout << "  Edge " << eid
            << " : q" << edge.getsourceid()
            << " -- q" << edge.gettargetid()
            << "  weight=" << edge.getweight()
            << "\n";
    }

    cout << "===================================================\n\n";
}
