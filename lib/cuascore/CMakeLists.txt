add_library(cuascore
  src/specialgradient.cpp
  src/CUASModel.cpp
  src/CUASArgs.cpp
  src/systemmatrix.cpp
  src/CUASSolver.cpp
  src/CUASFile.cpp
  src/SolutionHandler.cpp
  src/ModelReader.cpp
)

target_include_directories(cuascore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_property(TARGET cuascore PROPERTY CXX_STANDARD 14)

target_link_libraries(cuascore petscwrapper)

add_netcdf(cuascore)

add_cxxopts(cuascore)
add_spdlog(cuascore)

install(TARGETS cuascore DESTINATION lib)
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/ DESTINATION include FILES_MATCHING PATTERN "*.h")

