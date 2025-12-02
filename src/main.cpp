//main.cpp
//
/* ------This is main cpp for qcss ur ------- */

#include "circuit.h"
#include "QASMtoken.hpp"
#include "QASMscanner.hpp"
#include "QASMparser.h"
#include "mymeasure.h"

    //execute each procedure
#define parsing_ 1
#define gen_qpu  1
#define init_map 1
#define main_map 1
#define output_  1

using namespace std;
using namespace Qcircuit;

Qcircuit::QMapper mapper;
static CMeasure measure;

int main_2(int argc, char** argv){
    cout << endl;
    cout << "================================================================" << endl;
    cout << "         Multi-QPU Mapping by Ajou University QCSS Lab          " << endl;
    cout << "================================================================" << endl;

    measure.start_clock();

    #if parsing_
    mapper.parsing(argc, argv);
    measure.stop_clock("PARSING");
    #endif

    #if gen_qpu
    int num_qpu     = 3;
    int num_qubits  = 9;
    bool fully_connected_qubit  = true;
    bool fully_connected_qpu    = false;
    bool buffer_insertion       = true;
    mapper.generate_multi_qpu(num_qpu, num_qubits, fully_connected_qubit, fully_connected_qpu, buffer_insertion);
    measure.stop_clock("GENERATE MULTI-QPU");
    #endif

    #if init_map
    mapper.initial_mapping(num_qpu,num_qubits);
    measure.stop_clock("INITIAL MAPPING");
    #endif

    #if main_map
    mapper.main_mapping(mapper.Dgraph);
    measure.stop_clock("MAIN MAPPING");
    #endif

    #if output_
    mapper.FinalCircuit_info(mapper.FinalCircuit, true);
    mapper.write_output(mapper.FinalCircuit);
    measure.stop_clock("END PROGRAM");
    #endif
   
    cout << endl << "*** Time & Memory Measure ***" << endl;
    measure.print_clock();
    measure.printMemoryUsage();

    return 0;
}
