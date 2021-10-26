include(ExternalProject)

# Note, the spdlog library is header only.
if(DEFINED SPDLOG_INCLUDE)
  add_custom_target(spdlog
          COMMAND echo "SPDLOG_INCLUDE predefined: ${SPDLOG_INCLUDE}, skipping rebuild")
else()
  # in lib/petscwrapper/include/Logger.h: #include "spdlog/spdlog.h"
  find_path(SPDLOG_INCLUDE spdlog/spdlog.h)
  if(SPDLOG_INCLUDE)
    add_custom_target(spdlog
            COMMAND echo "SPDLOG_INCLUDE found in: ${SPDLOG_INCLUDE}, skipping rebuild")
  else()
    message("SPDLOG library not found, download into extern during make")
    ExternalProject_Add(spdlog
      SOURCE_DIR          ${CMAKE_CURRENT_SOURCE_DIR}/extern/spdlog
      GIT_REPOSITORY      "https://github.com/gabime/spdlog.git"
      GIT_TAG             v1.x
      CONFIGURE_COMMAND   ""
      BUILD_COMMAND       ""
      INSTALL_COMMAND     ""
      TEST_COMMAND        ""
      GIT_SHALLOW         true
    )
    set(SPDLOG_INCLUDE ${CMAKE_CURRENT_SOURCE_DIR}/extern/spdlog/include)
  endif()
endif()

function(add_spdlog target)
  add_dependencies(${target}
    spdlog
  )

  target_include_directories(${target} SYSTEM PUBLIC
    ${SPDLOG_INCLUDE}
  )
endfunction()

