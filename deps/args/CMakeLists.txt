set(ARGS_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
message("Using Args in ${ARGS_INCLUDE_DIR}")
add_library(args INTERFACE)
target_include_directories(args INTERFACE ${ARGS_INCLUDE_DIR})
