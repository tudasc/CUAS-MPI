# Doxygen + Swinx + Breathe

include(UseBreathe)

add_breathe_doc(
  SOURCE_DIR
  ${CMAKE_CURRENT_SOURCE_DIR}/
  BUILD_DIR
  ${CMAKE_CURRENT_BINARY_DIR}/_build
  CACHE_DIR
  ${CMAKE_CURRENT_BINARY_DIR}/_doctrees
  HTML_DIR
  ${CMAKE_CURRENT_BINARY_DIR}/html
  DOXY_FILE
  ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile.in
  CONF_FILE
  ${CMAKE_CURRENT_SOURCE_DIR}/conf.py.in
  TARGET_NAME
  docs
  COMMENT
  "HTML documentation")
