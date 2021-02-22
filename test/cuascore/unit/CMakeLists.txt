add_executable(test-nodata
        noData.cpp
        )
target_link_libraries(test-nodata cuascore)
target_link_libraries(test-nodata gtest)
set_property(TARGET test-nodata PROPERTY CXX_STANDARD 14)

add_executable(test-fillMatrix
  fill_matrix_test.cpp
)
target_link_libraries(test-fillMatrix cuascore)
target_link_libraries(test-fillMatrix gtest)
set_property(TARGET test-fillMatrix PROPERTY CXX_STANDARD 14)


add_executable(test-helper
  helperTest.cpp
)
target_link_libraries(test-helper cuascore)
target_link_libraries(test-helper gtest)
set_property(TARGET test-helper PROPERTY CXX_STANDARD 14)

add_executable(test-cxxargs
        cxxoptsTest.cpp
)
target_link_libraries(test-cxxargs cuascore)
target_link_libraries(test-cxxargs gtest)
set_property(TARGET test-cxxargs PROPERTY CXX_STANDARD 14)

add_executable(setup-test
  SetupTest.cpp
)
target_link_libraries(setup-test cuascore)
set_property(TARGET setup-test PROPERTY CXX_STANDARD 14)

# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

add_executable(test-CUASKernels
        CUASKernelsTest.cpp
        )
target_link_libraries(test-CUASKernels cuascore)
target_link_libraries(test-CUASKernels gtest)
set_property(TARGET test-CUASKernels PROPERTY CXX_STANDARD 14)

add_executable(test-gradient
        GradientTest.cpp
        )
target_link_libraries(test-gradient cuascore)
target_link_libraries(test-gradient gtest)
set_property(TARGET test-gradient PROPERTY CXX_STANDARD 14)
