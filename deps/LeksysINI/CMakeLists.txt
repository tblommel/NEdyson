set(INIPARSER_INCLUDE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)
message("Using INIParser in ${INIPARSER_INCLUDE_DIR}")
add_library(iniparser INTERFACE)
target_include_directories(iniparser INTERFACE ${INIPARSER_INCLUDE_DIR})
