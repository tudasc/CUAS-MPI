cmake_minimum_required(VERSION 3.3)

project(CUAS-MPI)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PETSc REQUIRED)

find_package(Git)
include(gitversion)
# Add a custom command that produces version.cpp, plus
# a dummy output that's not actually produced, in order
# to force gitversion.cmake to always be re-run before the build
ADD_CUSTOM_COMMAND(
        OUTPUT ${PROJECT_BINARY_DIR}/version.cpp
        ${PROJECT_BINARY_DIR}/_version.cpp
        COMMAND ${CMAKE_COMMAND} -P
        ${PROJECT_SOURCE_DIR}/cmake/gitversion.cmake)

find_package(NetCDF REQUIRED)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

include(ToolchainOptions)

add_subdirectory(lib)
add_subdirectory(tools)

option(BUILD_TESTS "enables building tests" OFF)

if(BUILD_TESTS)
  add_subdirectory(test)
endif()

