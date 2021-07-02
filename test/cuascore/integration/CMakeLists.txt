add_executable(test-noData
        noDataTest.cpp
        )
target_link_libraries(test-noData cuascore)
target_link_libraries(test-noData gtest)
set_property(TARGET test-noData PROPERTY CXX_STANDARD 14)

install(TARGETS test-noData DESTINATION test/cuascore/integration)

