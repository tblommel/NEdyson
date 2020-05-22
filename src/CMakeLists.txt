add_library(libNEdyson
  molHFsolver.cpp molHFsolver.h
  dyson.cpp dyson.h
  dyson_omp.cpp dyson_omp.h
  elementops.cpp elementops.h
  function.cpp function.h
  greens.cpp greens.h
  greens_tstp.cpp greens_tstp.h
  hubbchainops.cpp hubbchainops.h
  integrals.cpp integrals.h
  integration.cpp integration.h
  mat_utils.cpp mat_utils.h
  spectral.cpp spectral.h
  vie2.cpp vie2.h
  utils.cpp utils.h
  NEdyson.h
  param.h
)

target_include_directories(libNEdyson INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(libNEdyson PUBLIC 
  DATA_DIR="${PROJECT_SOURCE_DIR}/../NEdyson/data/")
target_link_libraries(libNEdyson PUBLIC
  alps::alps-gf
  Eigen3::Eigen
)
