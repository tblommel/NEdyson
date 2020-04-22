cmake_minimum_required(VERSION 3.1)

project(MySolver CXX)

find_package(Eigen3 REQUIRED)

message(STATUS ${EIGEN3_INCLUDE_DIR})

add_executable(main hubb_chain_2b.cpp headers2/integrals.cpp headers2/spectral.cpp headers2/greens.cpp headers2/vie2.cpp headers2/greens_tstp.cpp headers2/integration.cpp headers2/dyson.cpp headers2/function.cpp  headers2/elementops.cpp headers2/mat_utils.cpp)
target_include_directories(main PUBLIC ${EIGEN3_INCLUDE_DIR} .)
