cmake_minimum_required(VERSION 3.3)

project(CUAS-MPI)

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake)

find_package(PETSc REQUIRED)
message (STATUS "PETSC COMPILER ${PETSC_COMPILER}")
if (DEFINED PETSC_VERSION)
  message (STATUS "PETSC VERSION ${PETSC_VERSION}")
  if (PETSC_VERSION VERSION_LESS 3.9)
    # Force to look for PETSc again if the version we just found is too old:
    set(PETSC_CURRENT "OFF" CACHE BOOL "" FORCE)
    # Stop with an error message.
    message(FATAL_ERROR "CUAS-MPI requires PETSc version 3.9 or newer (found ${PETSC_VERSION}).")
  endif()
endif (DEFINED PETSC_VERSION)

# check for PETSc features MUMPS and PARMETIS in petscconf.h
MACRO(_petsc_feature_check _var _regex)
  FILE(STRINGS "${PETSC_INCLUDE_CONF}/petscconf.h" PETSC_${_var}_STRING
          REGEX "${_regex}")
  IF("${PETSC_${_var}_STRING}" STREQUAL "")
    SET(PETSC_WITH_${_var} FALSE)
  ELSE()
    SET(PETSC_WITH_${_var} TRUE)
    message (STATUS "PETSC_WITH_${_var} " TRUE)
  ENDIF()
ENDMACRO()
_petsc_feature_check(MUMPS "#define.*PETSC_HAVE_MUMPS 1")
_petsc_feature_check(PARMETIS "#define.*PETSC_HAVE_PARMETIS 1")

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

option(BUILD_TESTS "Build test programs" ON)

if(BUILD_TESTS)
  add_subdirectory(test)
endif()

