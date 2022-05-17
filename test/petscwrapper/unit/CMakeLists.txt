# Now simply link against gtest or gtest_main as needed. Eg
#add_executable(example example.cpp)
#target_link_libraries(example gtest_main)
#add_test(NAME example_test COMMAND example)

add_executable(test-PETScVector
        PETScVectorTest.cpp
        )
target_link_libraries(test-PETScVector petscwrapper)
target_link_libraries(test-PETScVector gtest)

add_executable(test-PETScMatrix
        PETScMatrixTest.cpp
        )
target_link_libraries(test-PETScMatrix cuascore)
target_link_libraries(test-PETScMatrix gtest)

add_executable(test-PETScSolver
        PETScSolverTest.cpp
        )
target_link_libraries(test-PETScSolver petscwrapper)
target_link_libraries(test-PETScSolver gtest)

add_executable(test-PETScGrid
        PETScGridTest.cpp
        )
target_link_libraries(test-PETScGrid petscwrapper)
target_link_libraries(test-PETScGrid gtest)

add_executable(test-fillgrid
        fillgridTest.cpp
        )
target_link_libraries(test-fillgrid petscwrapper)
target_link_libraries(test-fillgrid gtest)

install(TARGETS test-PETScVector DESTINATION test/petscwrapper/unit)
install(TARGETS test-PETScMatrix DESTINATION test/petscwrapper/unit)
install(TARGETS test-PETScSolver DESTINATION test/petscwrapper/unit)
install(TARGETS test-PETScGrid DESTINATION test/petscwrapper/unit)
install(TARGETS test-fillgrid DESTINATION test/petscwrapper/unit)
