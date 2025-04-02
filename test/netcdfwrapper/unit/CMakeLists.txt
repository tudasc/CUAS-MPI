# Now simply link against gtest or gtest_main as needed. Eg
# add_executable(example example.cpp) target_link_libraries(example gtest_main)
# add_test(NAME example_test COMMAND example)

add_executable(test-NetCDFFile NetCDFFileTest.cpp)
target_link_libraries(test-NetCDFFile netcdfwrapper)
target_link_libraries(test-NetCDFFile gtest)

install(TARGETS test-NetCDFFile DESTINATION test/netcdfwrapper/unit)
