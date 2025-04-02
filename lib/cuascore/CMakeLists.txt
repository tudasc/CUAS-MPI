add_library(
  cuascore
  src/CUASModel.cpp
  src/CUASArgs.cpp
  src/systemmatrix.cpp
  src/CUASSolver.cpp
  src/SolutionHandler.cpp
  src/ModelReader.cpp
  src/utilities.cpp
  src/version.cpp)

target_include_directories(cuascore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(cuascore petscwrapper)
target_link_libraries(cuascore netcdfwrapper)

add_dependencies(cuascore gitversion)

add_cxxopts(cuascore)
add_yamlcpp(cuascore)

add_spdlog(cuascore)

install(TARGETS cuascore DESTINATION lib)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.h")
