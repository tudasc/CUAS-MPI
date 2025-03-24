# derived from https://github.com/jbeder/yaml-cpp

include(FetchContent)

FetchContent_Declare(
  yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git
  GIT_TAG 0.8.0)

FetchContent_MakeAvailable(yaml-cpp)

if(TARGET yaml-cpp)
  message(STATUS "yaml-cpp added SUCCESS")
else()
  message(FATAL_ERROR "yaml-cpp not added")
endif()

get_target_property(YAMLCPP_INCLUDE yaml-cpp INTERFACE_INCLUDE_DIRECTORIES)
message(STATUS "YAMLCPP_INCLUDE ${YAMLCPP_INCLUDE}")

# the following environment variables are not found but yaml-cpp still works
# get_target_property(YAMLCPP_LIB_DIR yaml-cpp IMPORTED_LOCATION) message(STATUS
# "YAMLCPP_LIB_DIR ${YAMLCPP_LIB_DIR}") get_target_property(YAMLCPP_LIB yaml-cpp
# LOCATION) message(STATUS "YAMLCPP_LIB ${YAMLCPP_LIB}")

function(add_yamlcpp target)
  # add_dependencies(${target} yaml-cpp::yaml-cpp)

  # target_include_directories(${target} SYSTEM PUBLIC ${YAMLCPP_INCLUDE})

  target_link_libraries(${target} yaml-cpp::yaml-cpp)

  # we would like to use 'target_link_libraries(${target} PUBLIC
  # yaml-cpp::yaml-cpp)' but Keyword-Syntax and Plain-Syntax cannot be combined
  # in one project
endfunction()
