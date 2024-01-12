add_library(libNEdyson
  molHFSolver.cpp molHFSolver.h
  SCHFSolver.cpp SCHFSolver.h
  mol2bSolver.cpp mol2bSolver.h
  dyson.cpp dyson.h
  greens.cpp greens.h
  integration.cpp integration.h
  utils.cpp utils.h
  molNEsim.cpp molNEsim.h molNEsim_impl.h molNEsim_decomp_impl.h molNEsim_spin_impl.h molNEsim_decomp_spin_impl.h
  params.cpp params.h
)



target_include_directories(libNEdyson INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_definitions(libNEdyson PUBLIC 
  DATA_DIR="${DATA_DIR}"
  TIMING_DATA_DIR="${TIMING_DATA_DIR}")
target_link_libraries(libNEdyson PUBLIC
  alps::alps-gf
  Eigen3::Eigen
#  HighFive
  libgfmol
  cpsitop
  OpenMP::OpenMP_CXX
)
