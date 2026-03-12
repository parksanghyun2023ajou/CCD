cmake_minimum_required(VERSION 3.20)

set(CMAKE_CXX_STANDARD 17)

project(CCD VERSION 0.1)

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
    src/main.cpp
    src/run.cpp
    src/background_main_mapping.cpp
    src/mapping_preparation.cpp
    )

add_executable(CCD src/run.cpp ${SRC} ${TOOL})
target_include_directories(CCD PUBLIC ${CMAKE_SOURCE_DIR}/include )
target_link_libraries(CCD PRIVATE metis)
set_target_properties(
    CCD
    PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/tool
)
