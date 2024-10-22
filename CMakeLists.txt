cmake_minimum_required(VERSION 3.11)

project(NEdyson)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)
set(CMAKE_CXX_STANDARD 14)


find_package(ALPSCore REQUIRED COMPONENTS gf)
find_package(Boost REQUIRED)
message(STATUS ${Boost_INCLUDE_DIRS})
find_package(Eigen3 REQUIRED)
message(STATUS ${EIGEN3_INCLUDE_DIR})
find_package(HDF5 REQUIRED COMPONENTS C)
find_package(OpenMP REQUIRED)

add_subdirectory(${PROJECT_SOURCE_DIR}/deps/args)
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/LeksysINI)

#fetch HighFive from GitHub
include(CPM)
CPMAddPackage(
    NAME HighFive
    GITHUB_REPOSITORY BlueBrain/HighFive
    GIT_TAG master
    OPTIONS
    "HIGHFIVE_USE_BOOST OFF"
    "HIGHFIVE_USE_EIGEN ON"
    "HIGHFIVE_UNIT_TESTS OFF"
    "HIGHFIVE_EXAMPLES OFF"
)

include_directories(${INCLUDE_DIR})
link_directories(${LINK_DIR})

include_directories("${PROJECT_SOURCE_DIR}/deps/args/include/")
include_directories("${PROJECT_SOURCE_DIR}/deps/LeksysINI/include/")

add_subdirectory(src)
add_subdirectory(programs)

enable_testing()
add_subdirectory(${PROJECT_SOURCE_DIR}/deps/catch2)
add_subdirectory(tests)
