cmake_minimum_required(VERSION 3.12)

project(NEdyson)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 11)

option(omp "Build with OpenMP functionality" OFF)
if(omp)
  message(STATUS "Build with OpenMP functionality")
  find_package(OpenMP REQUIRED)
  add_definitions("-DUSE_OMP")
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

find_package(ALPSCore REQUIRED COMPONENTS gf)
find_package(Eigen3 REQUIRED)
message(STATUS ${EIGEN3_INCLUDE_DIR})
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIR})

include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/src/)

add_subdirectory(src)

add_subdirectory(programs)

add_subdirectory(CalcA)
