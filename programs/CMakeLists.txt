add_executable(molNEdyson main_mol.cpp)
target_include_directories(molNEdyson INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(molNEdyson PUBLIC libNEdyson
                                        libgfmol
                                        args
                                        iniparser)
