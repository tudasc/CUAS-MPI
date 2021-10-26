include(ExternalProject)

# Note, the cxxopts library is header only.
if(DEFINED CXXOPTS_INCLUDE)
  add_custom_target(cxxopts
          COMMAND echo "CXXOPTS_INCLUDE predefined: ${CXXOPTS_INCLUDE}, skipping rebuild")
else()
  # in lib/cuascore/src/CUASArgs.cpp: #include "cxxopts.hpp""
  find_path(CXXOPTS_INCLUDE cxxopts.hpp)
  if(CXXOPTS_INCLUDE)
    add_custom_target(cxxopts
            COMMAND echo "CXXOPTS_INCLUDE found in: ${CXXOPTS_INCLUDE}, skipping rebuild")
  else()
    message("CXXOPTS library not found, download into extern during make")
    ExternalProject_Add(cxxopts
      SOURCE_DIR          ${CMAKE_CURRENT_SOURCE_DIR}/extern/cxxopts
      GIT_REPOSITORY      "https://github.com/jarro2783/cxxopts.git"
      GIT_TAG             v2.2.1
      CONFIGURE_COMMAND   ""
      BUILD_COMMAND       ""
      INSTALL_COMMAND     ""
      TEST_COMMAND        ""
      GIT_SHALLOW         true
    )
    set(CXXOPTS_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/extern/cxxopts/include)
  endif()
endif()

function(add_cxxopts target)
  add_dependencies(${target}
    cxxopts
  )

  target_include_directories(${target} SYSTEM PUBLIC
    ${CXXOPTS_INCLUDE}
  )
endfunction()
