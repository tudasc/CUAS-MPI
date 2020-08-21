add_executable(vec-test
  VecTest.cpp
)

target_link_libraries(vec-test cuascore)
target_include_directories(vec-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)
