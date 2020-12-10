add_executable(vec-test
  VecTest.cpp
)
target_link_libraries(vec-test cuascore)

add_executable(mat-test
  MatTest.cpp
)
target_link_libraries(mat-test cuascore)

add_executable(grid-test
  GridTest.cpp
)
target_link_libraries(grid-test cuascore)



# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

add_executable(test-petscsolver
        SolverTest.cpp
        )
target_link_libraries(test-petscsolver petscwrapper)
target_link_libraries(test-petscsolver gtest)
set_property(TARGET test-petscsolver PROPERTY CXX_STANDARD 14)
