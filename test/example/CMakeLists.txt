# Now simply link against gtest or gtest_main as needed. Eg
add_executable(example example.cpp)
target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

add_executable(mpiexample mpiexample.cpp)
target_link_libraries(mpiexample gtest)
#add_test(NAME example_test COMMAND example)
