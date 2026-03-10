//circuit.h

#ifndef __CIRCUIT__
#define __CIRCUIT__ 0

//io
#include "metis.h"
#include <iostream>
#include <iomanip>
#include <fstream>
//ds
#include <vector>
#include <string>
#include <map>
#include <queue>
#include <list>
#include <set>
#include <unordered_map>
#include <unordered_set>
//cpp
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <cfloat>
#include <climits>
#include <cassert>
//op
#include <algorithm>
#include <iterator>
#include <utility>
#include <functional>
#include <random>

#define QASM_PI 3.141592653589793238463

using namespace std;

namespace Qcircuit
{
    enum class GATETYPE
    {
        U, RX, RZ, CNOT, X, Y, Z, H, S, SDG, T, TDG, MEASURE
    };

    //Circuit info
    struct Gate
    {
        int id;
        int control, target;
        GATETYPE type;
        char output_type[128];
    };

    struct Circuit
    {
        vector<Gate> nodeset;
    };
    
    struct Node;
    struct Edge;

    struct Node
    {
        private:
            ////variables
            static int global_id;               // 바뀌지 않는 글로벌 id 설정
            int id;                             // Local 인덱스 설정
            int qpu_id;                         // QPU id 설정
            int weight;                         // 가중치 설정

        public:
            set<int> edges;                     // 정수를 저장하는 set

        public:
            ////constructors
            Node();                             // 인자 없는 Node 함수 호출
            Node(int weight);                   // 인자가 가중치인 Node 함수

            ////getters
            int getglobalid();                  // 해당 노드에 글로벌 아이디를 생성하는 함수 선언
            int getid();                        // id를 얻는 함수 선언
            int getweight();                    // weight를 얻는 함수 선언
            void setweight(int weight);         // weight를 설정하는 함수 선언

            static void init_global_id();       // 구조체의 객체를 설정하지 않고도 호출할 수 있음, 정적 멤버에만 접근 가능
    };
    struct Edge
    {
        private:
            ////variables
            static int global_id;
            int id;
            double weight;
            public:
            Node* source;                       // Node형 pointer source 선언
            Node* target;                       // Node형 pointer target 선언

        public:
            ////constructors
            Edge();                            // 인자 없는 Edge 함수 호출
            Edge(Node* source, Node* target, double weight);
            // 연결된 node 정보와 weight와 함께 edge 함수 선언

            ////getters
            int getglobalid() const;
            int getid()const;
            double getweight()const;
            int getsourceid()const;                  // source의 id 얻는 함수 선언
            int gettargetid()const;                  // target의 id 얻는 함수 선언
            void setweight(double weight);
            void setsource(Node* source);       // source를 얻는 함수 선언
            void settarget(Node* target);       // target을 얻는 함수 선언

            static void init_global_id();

    };

    struct Graph
    {
        ////graph node and edge
        map<int, Node> nodeset;                 // node를 저장하는 map
        map<int, Edge> edgeset;                 // edge를 저장하는 map
        
        map<int, int> nodeid;       //table index <--> graph node id
        
        int node_size;                          // node의 사이즈 저장 변순
        int** dist;                             // 이중 포인터 선언
        int center;
        vector<pair<int,int>> degree_info;
        queue<pair<int,int>> degree_queue;                             // 센터 선언

        public:
            ////constructors
            Graph();

            ////methods
            //graph.cpp
            void init_global_id();
            int  addnode(int weight = 1);
            bool deletenode(int id);
            bool addedge(int source, int target, double weight);
            bool deleteedge(int id);

            //graphFunction.cpp
            void build_dist_table();
            void print_dist_table();
            void delete_dist_table();
            int bfs(int start, int end);
            int bfs_weight(int start, int end);
            int cal_graph_center();
            pair<int, bool> graph_characteristics(bool i);
    };

    //for mapping
    class QMapper
    {
        public:
        ///////////// variables
            //file names
            string fileName_input;              // 문자열 파일 이름 input 변수 선언
            string fileName_output;             // 문자열 파일 이름 output 변수 선언

            //number of qubits;
            unsigned int nqubits;
            unsigned int positions;
            Graph multi_qpu_graph;
            Graph qpu_degree_graph;
            Graph interactionGraph;

            //Quantum circuit
            Circuit Dgraph;                     // Gate 타입의 벡터 nodeset 객체
            Circuit Dgraph_cnot;
            Circuit FinalCircuit;
              //Dlist
            vector<list<int> > Dlist;
            vector<list<int> > Dlist_all;

               // Gate 타입의 벡터 nodeset 객체
            int num_qubits;
             int node_id;
            int add_cnot_num;
            int add_swap_num;
            int add_bridge_num;
            int least_cnot_num;     //least cnot_num
            int least_swap_num;     //least swap_num
            int least_bridge_num;   //least bridge_num
            double param_alpha;
            int param_beta;
            vector<vector<double>> cross_table;
            vector<idx_t> partition_result;
            vector<Graph> Subsets;
            vector<pair<int,pair<int,Edge>>> inter_Sub_edge;
            vector <pair<Graph,int>> matching_info;
            vector<map<int, int>> layout_L;  //L[logical_qubit] = physical_qubit
            vector<map<int, int>> qubit_Q;   //Q[physical_qubit] = logical_qubit
            //output
            int add_2q_num;
            float fidelity;

        //////////// functions
            //parser.cpp
            void parsing(int argc, char** argv);
            void QASM_parse(string filename);

            //generate_qpu.cpp
            void generate_multi_qpu(int num_qpu, int num_qubits, 
                bool fully_connected_qubit, bool fully_connected_qpu, bool buffer_insertion);
                
           //interaction_graph.cpp
           void make_CNOT(bool i);
           Graph make_interactionGraph(bool i);
           Graph make_interactionNumberGraph(bool i);
           Graph make_interactionMixgraph(bool i, int n);
           void print_interactionGraph( Graph& g);
           int extract_qpu_idx(int logical_qubit);
           //METIS_executor.cpp
           void run_metis_partition(Graph& interactionGraph, int num_parts);
           void analyze_cross_partition_edges(Graph& g, const vector<idx_t>& part, int num_qpu);
           void reconstruct_info(Graph& interactiongraph, int num_parts);
           void matching();
           // cicuit_process.cpp
           void circuit_processing(int num_qpu);


                //initial_mapping.cpp
            int physical_degree_local(int qpu, int local_idx, int num_qubits);
            vector<int> extract_logicals_from_subgraph(const Graph& g);
            void initial_mapping(int num_qpu,int num_qubit);
            // Mapping_Function.cpp
            int cal_SWAP_effect(const int q1, const int q2, const int Q1, const int Q2);
            double cal_MCPE(const pair<int, int> p, Circuit& dgraph);
            void find_max_cost(pair<int, int>& SWAP, int& gateid, double& max_cost,
                                             vector< pair< pair<int, int>, pair<int, double> > >& v, list<int>& act_list,
                                             vector< pair< pair<int, int>, pair<int, double> > >& history);
            void layout_swap(const int b1, const int b2);
            void add_swap(int b1, int b2, Circuit& graph);
            void add_cnot(int c_qubit, int t_qubit, Circuit& graph);
            void add_bridge(int qs, int qt, Circuit& graph);
            int sort_degree_return(vector<int>& candi_loc, Graph& graph,int num_qpu);
            void sort_degree(vector<int>& candi_loc, Graph& graph,int num_qpu);
            void make_Dlist(Circuit& dgraph);
            void make_Dlist_all(Circuit& dgraph);
            void bfs_queue_gen(queue<int>& queue, Graph& graph, int start, bool order);
            void make_ref_loc(vector<int>& ref_loc, Graph& graph, int start, bool order);
            void make_candi_loc(int current_q, Graph& graph, vector<int>& candi_loc_1, vector<int>& candi_loc_2, int& degree,int num_qpu);
            void make_candi_loc_dist(int num_qpu,int current_qc, vector<int>& candi_loc, vector<int>& ref_loc, Graph& coupling_graph, Graph& interaction_graph, bool equal_order);
            void update_front_n_act_list(list<int>& front_list, list<int>& act_list, vector<bool>& frozen);
            void find_singlequbit_list(int gateid, list<int>& singlequbit_list, Circuit& dgraph);
            void update_act_dist2_list(list<int>& act_dist2_list, list<int>& act_list, Circuit& dgraph);
            bool check_direct_act_list(list<int>& act_list, list<int>& singlequbit_list, vector<bool>& frozen, Circuit& dgraph);
            void generate_candi_list(list<int>& act_list, vector< pair< pair<int, int>, int> >& candi_list, Circuit& dgraph);


            //main_mapping.cpp
            void main_mapping(Circuit& dgraph);
            void mapping_machine(bool cost_flag, const pair<int, int> c,  Circuit& dgraph, const int q1, const int q2, const int Q1, const int Q2);

            //outputwriter.cpp
            void FinalCircuit_info(Circuit& graph, bool final_circuit);
            void write_output(Circuit& graph);
    };
}

using namespace Qcircuit;
    //// Circuit
inline ostream& operator << (ostream& os, const GATETYPE& type)
{
    switch(type)
    {
        case GATETYPE::U    : return os << "u"  ;
        case GATETYPE::RX   : return os << "rx" ;
        case GATETYPE::RZ   : return os << "rz" ;
        case GATETYPE::CNOT : return os << "cx" ;
        case GATETYPE::X    : return os << "x"  ;
        case GATETYPE::Y    : return os << "y"  ;
        case GATETYPE::Z    : return os << "z"  ;
        case GATETYPE::H    : return os << "h"  ;
        case GATETYPE::S    : return os << "s"  ;
        case GATETYPE::SDG  : return os << "sdg";
        case GATETYPE::T    : return os << "t"  ;
        case GATETYPE::TDG  : return os << "tdg";
        case GATETYPE::MEASURE  : return os << "measure";
        default             : return os;
    }
}

inline ostream& operator << (ostream& os, const Gate& gate)
{
    os << "Gate";
    if(gate.control == -1) //U, RX, RZ gate
        os << "(" << setw(5) << gate.id << ")\t" << setw(2) << gate.type << " " << "q[" << gate.target << "]";
    else //CNOT gate
        os << "(" << setw(5) << gate.id << ")\t" << setw(2) << gate.type << " " << "q[" << gate.control << "], q[" << gate.target << "]";
    return os;   
}

inline ostream& operator << (ostream& os, const vector<Gate>& nodeset)
{
    for(auto gate : nodeset)
        os << gate << endl;
    return os;
}

inline ostream& operator << (ostream& os, const list<int>& list)
{
    for(auto v : list)
        os << v << " ";
    return os;
}

#endif 