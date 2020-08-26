cmake_minimum_required(VERSION 3.1)

add_executable(calcA.x calcA.cpp)
target_link_libraries(calcA.x libNEdyson)
add_executable(calcA_tti.x calcA_tti.cpp)
target_link_libraries(calcA_tti.x libNEdyson)
