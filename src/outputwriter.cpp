//outputwriter.cpp
#include <iomanip>
#include "circuit.h"

#define PRINT_OUTPUTWRITER 0

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::FinalCircuit_info(Circuit& graph, bool final_circuit)
{
    if(final_circuit)
    {
        cout << "*** Final Circuit ***" << endl;
        cout << "# add_2q_num: "    << add_2q_num << endl;
        cout << "# fidelity: " << fixed << setprecision(6) << fidelity << endl;
    }
    #if PRINT_OUTPUTWRITER
    //cout << graph.nodeset << endl;
    cout << " ==> Print Final Circuit Successfully Finished" << endl;
    #endif
}

void Qcircuit::QMapper::write_output(Circuit& graph)
{
    ofstream of(fileName_output);
    of << "OPENQASM 2.0;" << endl;
    of << "include \"qelib1.inc\";" << endl;
    of << "qreg q[" << nqubits << "]" << endl;
    of << "creg c[" << nqubits << "]" << endl;
    for( auto g : graph.nodeset)
    {
        if(g.control==-1)
            of << g.type <<  g.output_type << " " << "q[" << g.target << "];" << endl;
        else
            of << g.type << " " << "q[" << g.control << "],q[" << g.target << "];" << endl; 
    }
    of.close();
    #if PRINT_OUTPUTWRITER
    cout << " ==> Output File (.qasm) Is Successfully Generated" << endl;
    #endif
}