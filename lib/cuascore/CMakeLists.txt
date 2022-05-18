add_library(cuascore
  src/specialgradient.cpp
  src/CUASModel.cpp
  src/CUASArgs.cpp
  src/systemmatrix.cpp
  src/CUASSolver.cpp
  src/NetCDFFile.cpp
  src/SolutionHandler.cpp
  src/ModelReader.cpp
  src/utilities.cpp
  src/version.cpp
)

target_include_directories(cuascore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(cuascore petscwrapper)

add_dependencies(cuascore gitversion)

add_netcdf(cuascore)

add_cxxopts(cuascore)
add_spdlog(cuascore)

install(TARGETS cuascore DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

