cmake_minimum_required(VERSION 3.3)

project(CUAS-MPI)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PETSc REQUIRED)
message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")

find_package(NetCDF REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(ToolchainOptions)

add_subdirectory(lib)
add_subdirectory(tools)

option(BUILD_TESTS "Build test programs" ON)

if(BUILD_TESTS)
  add_subdirectory(test)
endif()

