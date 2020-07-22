add_executable(runtests
  tests-main.cpp
  integrals-test.cpp
  tests-def.h
  dyson-test.cpp
  spectral-test.cpp
  free-test.cpp
  integration-test.cpp
  Nsim-test.cpp
  NTTIsim-test.cpp
  Dsim-test.cpp
  DTTIsim-test.cpp
  Ssim-test.cpp
  STTIsim-test.cpp
  SDsim-test.cpp
  SDTTIsim-test.cpp
  molHFSolver-test.cpp
  mol2bSolver-test.cpp)

link_directories("~/Libraries/gfmol/build/src/")
target_compile_definitions(runtests PRIVATE
  TEST_DATA_DIR="${PROJECT_SOURCE_DIR}/tests/data")
target_link_libraries(runtests PRIVATE
  catch2
  libNEdyson)

include(CTest)
include(Catch)

add_test(NAME NEdyson_tests
  COMMAND runtests -o report.xml -r junit)
