add_executable(test-noData noDataTest.cpp)
target_link_libraries(test-noData cuascore)
target_link_libraries(test-noData gtest)

install(TARGETS test-noData DESTINATION test/cuascore/integration)

add_executable(test-exactTransientSolution exactTransientSolutionTest.cpp)
target_link_libraries(test-exactTransientSolution cuascore)
target_link_libraries(test-exactTransientSolution gtest)

install(TARGETS test-exactTransientSolution
        DESTINATION test/cuascore/integration)
