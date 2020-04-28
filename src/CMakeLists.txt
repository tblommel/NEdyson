add_library(libNEdyson
  dyson.cpp dyson.h
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
