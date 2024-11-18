if(CUAS_ENABLE_DOCS STREQUAL "Doxygen")
  message(STATUS "Documentation mode: Doxygen")
  add_subdirectory(Doxygen)
else()
  message(
    WARNING "No valid documentation mode defined, using: Doxygen (fallback)")
  add_subdirectory(Doxygen)
endif()
