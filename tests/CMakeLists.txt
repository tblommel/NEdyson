add_executable(runtests
  tests-main.cpp
  tests-def.h
  dyson-test.cpp
  free-test.cpp
  integration-test.cpp
  hubb-test.cpp
  molHFSolver-test.cpp)

target_compile_definitions(runtests PRIVATE
  TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")
target_link_libraries(runtests PRIVATE
  catch2
  libNEdyson)

include(CTest)
include(Catch)

catch_discover_tests(runtests)
