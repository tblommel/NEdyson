cmake_minimum_required(VERSION 3.1)

project(CalcA)

find_package(Eigen3 REQUIRED)

message(STATUS ${EIGEN3_INCLUDE_DIR})

add_executable(main calcA.cpp ../headers2/spectral.cpp ../headers2/greens.cpp ../headers2/elementops.cpp)
target_include_directories(main PUBLIC ${EIGEN3_INCLUDE_DIR} ./..)
