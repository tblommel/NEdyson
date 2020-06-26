# List off executables to build
set( project_EXECUTABLES
     hubb_chain_2b.x
     hubb_chain_2b_omp.x
     hubb_chain_2b_B.x
     hubb_chain_tpp.x
)

# Set their sources
set( EXE_hubb_chain_2b.x_SOURCES hubb_chain_2b.cpp )
set( EXE_hubb_chain_2b_omp.x_SOURCES hubb_chain_2b_omp.cpp )
set( EXE_hubb_chain_2b_B.x_SOURCES hubb_chain_2b_B.cpp )
set( EXE_hubb_chain_tpp.x_SOURCES hubb_chain_tpp.cpp )

#foreach( exe ${project_EXECUTABLES} )
#    # declare the target...
#    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
#    target_link_libraries(${exe} libNEdyson
#                                 libgfmol)
#endforeach()

add_executable(molNEdyson main.cpp)
target_include_directories(molNEdyson INTERFACE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(molNEdyson PUBLIC libNEdyson
                                        libgfmol
                                        args
                                        iniparser)
