set(CATCH_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
message("Using Catch2 in ${CATCH_INCLUDE_DIR}")
add_library(catch2 INTERFACE)
target_include_directories(catch2 INTERFACE ${CATCH_INCLUDE_DIR})