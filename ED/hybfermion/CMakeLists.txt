cmake_minimum_required (VERSION 2.8.12)

project(alps_cthyb)
list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake)

include_directories(
    ${EIGEN3_INCLUDE_DIR}
)

# Request the ALPSCore package:
find_package(Eigen3 REQUIRED)

include_directories(${EIGEN3_INCLUDE_DIR})
add_executable(fermion hybfermion.cpp)
