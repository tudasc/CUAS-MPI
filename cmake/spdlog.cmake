include(ExternalProject)

if(DEFINED SPDLOG_INCLUDE)
  message("SPDLOG_INCLUDE predefined: ${SPDLOG_INCLUDE}")
else()
  find_path(SPDLOG_LIBRARY NAMES spdlog)
  if(SPDLOG_LIBRARY)
    set(SPDLOG_INCLUDE ${SPDLOG_LIBRARY}/spdlog/include)
    message("SPDLOG found in ${SPDLOG_INCLUDE}")
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

