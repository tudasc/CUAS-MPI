if(CUAS_ENABLE_DOCS STREQUAL "Doxygen")
  message(STATUS "Documentation mode: Doxygen")
  add_subdirectory(Doxygen)
elseif(CUAS_ENABLE_DOCS STREQUAL "DoxygenSphinxBreathe")
  message(STATUS "Documentation mode: DoxygenSphinxBreathe")
  add_subdirectory(DoxygenSphinxBreathe)
else()
  message(
    WARNING "No valid documentation mode defined, using: Doxygen (fallback)")
  add_subdirectory(Doxygen)
endif()
