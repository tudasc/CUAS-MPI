cmake_minimum_required(VERSION 3.16)

project(
  CUAS-MPI
  VERSION 0.2.0
  LANGUAGES CXX)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PETSc REQUIRED)
find_package(NetCDF REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

include(ToolchainOptions)

add_subdirectory(lib)
add_subdirectory(tools)

# The default options cache
option(CUAS_ENABLE_TESTS "Enables targets building tests." OFF)
option(CUAS_ENABLE_DOCS "Enables targets building documentation." OFF)

if(CUAS_ENABLE_TESTS)
  message(STATUS "creating targets to build tests is enabled")
  add_subdirectory(test)
endif()

if(CUAS_ENABLE_DOCS)
  message(STATUS "creating targets to build documentation is enabled")
  add_subdirectory(docs)
endif()
