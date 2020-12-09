add_executable(gradient-test
  GradientTest.cpp
)
target_link_libraries(gradient-test cuascore)

add_executable(nodata-test
  noData.cpp
)
target_link_libraries(nodata-test cuascore)

add_executable(fillMatrix-test
  fill_matrix_test.cpp
)
target_link_libraries(fillMatrix-test cuascore)

add_executable(helper-test
  helperTest.cpp
)
target_link_libraries(helper-test cuascore)

add_executable(args-test
  cxxoptsTest.cpp
)
target_link_libraries(args-test cuascore)



# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)
