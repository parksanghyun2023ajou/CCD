# include "circuit.h"

void Qcircuit::QMapper::circuit_processing(int num_qpu)
{
    make_CNOT(true);
    //Make interaction graph
    int CNOTNodesetSize = Dgraph_cnot.nodeset.size();
    int NodesetSize     = Dgraph.nodeset.size();
    param_beta=3;
    int n = CNOTNodesetSize/param_beta;
    if (nqubits<=10 || nqubits==20) 
        n = 1;
    if ( (nqubits==15 && NodesetSize<=5000) || nqubits==11)
        n = param_beta;
    interactionGraph = make_interactionMixgraph(true, n);
    interactionGraph.build_dist_table();

    
   // Graph interaction_graph_order     = make_interactionGraph(true);
    Graph interaction_graph_number    = make_interactionNumberGraph(true);
    //interaction_graph_order.build_dist_table();
    interaction_graph_number.build_dist_table();
    print_interactionGraph(interactionGraph);
    //run_metis_partition(interaction_graph_number,num_qpu);
    //reconstruct_info(interaction_graph_number,num_qpu);
    //analyze_cross_partition_edges(interaction_graph_number, partition_result, num_qpu);
    run_metis_partition(interactionGraph,num_qpu);
    reconstruct_info(interactionGraph,num_qpu);
    analyze_cross_partition_edges(interactionGraph,partition_result,num_qpu);
    matching();
}