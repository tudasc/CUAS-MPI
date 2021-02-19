add_library(cuascore
  src/specialgradient.cpp
  src/CUASModel.cpp
  src/CUASArgs.cpp
  src/fill_matrix_coo.cpp
  src/setup.cpp
  src/CUASSolver.cpp
)

target_include_directories(cuascore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)
set_property(TARGET cuascore PROPERTY CXX_STANDARD 14)

target_link_libraries(cuascore petscwrapper)

add_cxxopts(cuascore)

