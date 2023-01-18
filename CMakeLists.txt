cmake_minimum_required(VERSION 3.3)

project(CUAS-MPI)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PETSc REQUIRED)
find_package(NetCDF REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(ToolchainOptions)

add_subdirectory(lib)
add_subdirectory(tools)

option(BUILD_TESTS "enables building tests" OFF)

if(BUILD_TESTS)
  add_subdirectory(test)
endif()
