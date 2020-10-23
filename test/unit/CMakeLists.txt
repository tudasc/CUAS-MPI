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

add_executable(grid-test
  GridTest.cpp
)

target_link_libraries(grid-test cuascore)
target_include_directories(grid-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)

add_executable(gradient-test
  GradientTest.cpp
)

target_link_libraries(gradient-test cuascore)
target_include_directories(gradient-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)

add_executable(nodata-test
  noData.cpp
)

target_link_libraries(nodata-test cuascore)
target_include_directories(nodata-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)

add_executable(fillMatrix-test
  fill_matrix_test.cpp
)

target_link_libraries(fillMatrix-test cuascore)
target_include_directories(fillMatrix-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)

add_executable(helper-test
  helperTest.cpp
)

target_link_libraries(helper-test cuascore)
target_include_directories(helper-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)

add_executable(args-test
  cxxoptsTest.cpp
)

target_link_libraries(args-test cuascore)
target_include_directories(args-test PRIVATE ${PROJECT_SOURCE_DIR}/lib/include)
