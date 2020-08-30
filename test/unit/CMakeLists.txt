add_executable(vec-test
  VecTest.cpp
)

target_link_libraries(vec-test cuascore)
target_include_directories(vec-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)

add_executable(mat-test
  MatTest.cpp
)

target_link_libraries(mat-test cuascore)
target_include_directories(mat-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)
