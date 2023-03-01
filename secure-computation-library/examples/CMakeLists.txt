cmake_minimum_required( VERSION 3.14 )

project( scl_examples VERSION 0.3 DESCRIPTION "SCL example programs" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -O2 -std=gnu++17" )

include_directories( ../include )

set( libscl "${CMAKE_SOURCE_DIR}/libscl.so.0.3" )

## finite fields
set( finite_fields "01_finite_fields" )
add_executable( "${finite_fields}" "${finite_fields}.cc" )
target_link_libraries( "${finite_fields}" "${libscl}" pthread )

## matrices and vectors
set( matrix_and_vector "02_matrix_and_vector" )
add_executable( "${matrix_and_vector}" "${matrix_and_vector}.cc" )
target_link_libraries( "${matrix_and_vector}" "${libscl}" pthread )

## secret sharing
set( secret_sharing "03_secret_sharing" )
add_executable( "${secret_sharing}" "${secret_sharing}.cc" )
target_link_libraries( "${secret_sharing}" "${libscl}" pthread )

## simple networking
set( networking "04_networking" )
add_executable( "${networking}" "${networking}.cc" )
target_link_libraries( "${networking}" "${libscl}" pthread )

## network discovery
set( discovery "05_discovery" )
add_executable( "${discovery}" "${discovery}.cc" )
target_link_libraries( "${discovery}" "${libscl}" pthread )
