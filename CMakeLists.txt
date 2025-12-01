cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)

project(Convergent_Capstone_Design VERSION 0.1)

find_package(Boost REQUIRED)

set(TOOL
    src/qasm-tools/QASMtoken.cpp
    src/qasm-tools/QASMparser.cpp
    src/qasm-tools/QASMscanner.cpp)

set(SRC
    src/parser.cpp
    src/generate_qpu.cpp
    src/graph_function.cpp
    src/initial_mapping.cpp
    src/main_mapping.cpp
    src/outputwriter.cpp
    src/mymeasure.cpp
    src/graph.cpp
    src/interaction_graph.cpp
    src/metis_execution.cpp
    src/circuit_process.cpp
    )

add_executable(Convergent_Capstone_Design src/main.cpp ${SRC} ${TOOL})
target_include_directories(Convergent_Capstone_Design PUBLIC ${CMAKE_SOURCE_DIR}/include )
target_link_libraries(Convergent_Capstone_Design PRIVATE metis)
set_target_properties(
    Convergent_Capstone_Design
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tool
)
