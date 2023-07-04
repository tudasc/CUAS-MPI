add_executable(test-noData noDataTest.cpp)
target_link_libraries(test-noData cuascore)
target_link_libraries(test-noData gtest)

install(TARGETS test-noData DESTINATION test/cuascore/integration)

add_executable(test-exactTransientSolution exactTransientSolutionTest.cpp)
target_link_libraries(test-exactTransientSolution cuascore)
target_link_libraries(test-exactTransientSolution gtest)

install(TARGETS test-exactTransientSolution
        DESTINATION test/cuascore/integration)

add_executable(test-exactSteadySolution exactSteadySolutionTest.cpp)
target_link_libraries(test-exactSteadySolution cuascore)
target_link_libraries(test-exactSteadySolution gtest)

install(TARGETS test-exactSteadySolution DESTINATION test/cuascore/integration)

add_executable(test-exactCircularSteadySolution
               exactCircularSteadySolutionTest.cpp)
target_link_libraries(test-exactCircularSteadySolution cuascore)
target_link_libraries(test-exactCircularSteadySolution gtest)

install(TARGETS test-exactCircularSteadySolution
        DESTINATION test/cuascore/integration)
