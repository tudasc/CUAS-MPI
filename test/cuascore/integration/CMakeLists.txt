add_executable(test-noData
        noDataTest.cpp
        )
target_link_libraries(test-noData cuascore)
target_link_libraries(test-noData gtest)

install(TARGETS test-noData DESTINATION test/cuascore/integration)

