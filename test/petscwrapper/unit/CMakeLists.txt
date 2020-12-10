# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

add_executable(test-petscvec
        VecTest.cpp
        )
target_link_libraries(test-petscvec petscwrapper)
target_link_libraries(test-petscvec gtest)
set_property(TARGET test-petscvec PROPERTY CXX_STANDARD 14)

add_executable(test-petscmat
        MatTest.cpp
        )
target_link_libraries(test-petscmat cuascore)
target_link_libraries(test-petscmat gtest)
set_property(TARGET test-petscmat PROPERTY CXX_STANDARD 14)

add_executable(test-petscsolver
        SolverTest.cpp
        )
target_link_libraries(test-petscsolver petscwrapper)
target_link_libraries(test-petscsolver gtest)
set_property(TARGET test-petscsolver PROPERTY CXX_STANDARD 14)

add_executable(test-petscgrid
        GridTest.cpp
        )
target_link_libraries(test-petscgrid petscwrapper)
target_link_libraries(test-petscgrid gtest)
set_property(TARGET test-petscgrid PROPERTY CXX_STANDARD 14)
