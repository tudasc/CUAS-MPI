# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

add_executable(test-PETScVec
        PETScVecTest.cpp
        )
target_link_libraries(test-PETScVec petscwrapper)
target_link_libraries(test-PETScVec gtest)
set_property(TARGET test-PETScVec PROPERTY CXX_STANDARD 14)

add_executable(test-PETScMat
        PETScMatTest.cpp
        )
target_link_libraries(test-PETScMat cuascore)
target_link_libraries(test-PETScMat gtest)
set_property(TARGET test-PETScMat PROPERTY CXX_STANDARD 14)

add_executable(test-PETScSolver
        PETScSolverTest.cpp
        )
target_link_libraries(test-PETScSolver petscwrapper)
target_link_libraries(test-PETScSolver gtest)
set_property(TARGET test-PETScSolver PROPERTY CXX_STANDARD 14)

add_executable(test-PETScGrid
        PETScGridTest.cpp
        )
target_link_libraries(test-PETScGrid petscwrapper)
target_link_libraries(test-PETScGrid gtest)
set_property(TARGET test-PETScGrid PROPERTY CXX_STANDARD 14)
