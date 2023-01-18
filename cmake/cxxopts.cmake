include(ExternalProject)
# Note, the cxxopts library is header only.

if(EXISTS ${CXXOPTS_INCLUDE})
  message(STATUS "CXXOPTS_INCLUDE predefined, skipping rebuild")
  add_custom_target(cxxopts)
else()
  # in lib/cuascore/src/CUASArgs.cpp: #include "cxxopts.hpp""
  find_path(CXXOPTS_INCLUDE cxxopts.hpp
            HINTS ${CMAKE_CURRENT_SOURCE_DIR}/extern/cxxopts/include)
  if(EXISTS ${CXXOPTS_INCLUDE})
    message(STATUS "CXXOPTS_INCLUDE found, skipping rebuild")
    add_custom_target(cxxopts)
  else()
    message(
      STATUS "CXXOPTS_INCLUDE not found, download into extern during make")
    ExternalProject_Add(
      cxxopts
      SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/extern/cxxopts
      GIT_REPOSITORY "https://github.com/jarro2783/cxxopts.git"
      GIT_TAG v2.2.1
      CONFIGURE_COMMAND ""
      BUILD_COMMAND ""
      INSTALL_COMMAND ""
      TEST_COMMAND ""
      GIT_SHALLOW true)
    set(CXXOPTS_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/extern/cxxopts/include)
  endif()
endif()

message(STATUS "CXXOPTS_INCLUDE ${CXXOPTS_INCLUDE}")

function(add_cxxopts target)
  add_dependencies(${target} cxxopts)

  target_include_directories(${target} SYSTEM PUBLIC ${CXXOPTS_INCLUDE})
endfunction()
