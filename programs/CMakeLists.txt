# List off executables to build
set( project_EXECUTABLES
     hubb_chain_2b.x
     hubb_chain_tpp.x
     math_compare.x
)

# Set their sources
set( EXE_hubb_chain_2b.x_SOURCES hubb_chain_2b.cpp )
set( EXE_hubb_chain_tpp.x_SOURCES hubb_chain_tpp.cpp )
set( EXE_math_compare.x_SOURCES math_compare.cpp )

# Paths and Subdirs
include_directories(${CMAKE_INCLUDE_PATH})
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# Libraries
link_directories(${CMAKE_LIBRARY_PATH})

foreach( exe ${project_EXECUTABLES} )
    # declare the target...
    add_executable( ${exe} ${EXE_${exe}_SOURCES} )
    target_link_libraries(${exe} libNEdyson)
endforeach()
