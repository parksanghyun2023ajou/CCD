//parser.cpp
#include "circuit.h"
#include "QASMtoken.hpp"
#include "QASMscanner.hpp"
#include "QASMparser.h"

#define PRINT_PARSER 1
#define PRINT_PARSER_NODESET 1

using namespace std;
using namespace Qcircuit;

void Qcircuit::QMapper::parsing(int argc, char** argv)
{
    if(argc != 3)
    {
        cout << "ERROR: invalid arguments" << endl;
        exit(1);
    }

    fileName_input  = argv[1];
    fileName_output = argv[2];
    QASM_parse(fileName_input);
    
    #if PRINT_PARSER
    cout << "*** Option ***" << endl;
    cout << "Input file  : " << fileName_input  << endl;
    cout << "Output file : " << fileName_output << endl;
    cout << endl;
    #endif
}

void Qcircuit::QMapper::QASM_parse(string filename)
{
    Circuit* circuit;

    //Parsing input file
    circuit = &Dgraph;

    //Start Parser
    QASMparser* parser = new QASMparser(filename);
    parser->Parse();    // QASM parse
    nqubits = parser->getNqubits();
    circuit->nodeset = parser->getGatelists();
    #if PRINT_PARSER_NODESET
    cout << circuit->nodeset << endl;
    #endif

    delete parser;
}
