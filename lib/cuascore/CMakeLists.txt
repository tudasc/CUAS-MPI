add_library(cuascore
  src/specialgradient.cpp
  src/CUASModel.cpp
  src/parseCxxopts.cpp
)

target_include_directories(cuascore PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(cuascore petscwrapper)

add_cxxopts(cuascore)

