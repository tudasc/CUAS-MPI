add_executable(test-noData
        noDataTest.cpp
        )
target_link_libraries(test-noData cuascore)
target_link_libraries(test-noData gtest)
set_property(TARGET test-noData PROPERTY CXX_STANDARD 14)

add_executable(test-fillMatrix
        fillMatrixTest.cpp
)
target_link_libraries(test-fillMatrix cuascore)
target_link_libraries(test-fillMatrix gtest)
set_property(TARGET test-fillMatrix PROPERTY CXX_STANDARD 14)

add_executable(test-CUASKernels
        CUASKernelsTest.cpp
        )
target_link_libraries(test-CUASKernels cuascore)
target_link_libraries(test-CUASKernels gtest)
set_property(TARGET test-CUASKernels PROPERTY CXX_STANDARD 14)

add_executable(test-specialgradient
        specialgradientTest.cpp
        )
target_link_libraries(test-specialgradient cuascore)
target_link_libraries(test-specialgradient gtest)
set_property(TARGET test-specialgradient PROPERTY CXX_STANDARD 14)

add_executable(test-fillNoData
        fillNoDataTest.cpp
        )
target_link_libraries(test-fillNoData cuascore)
target_link_libraries(test-fillNoData gtest)
set_property(TARGET test-fillNoData PROPERTY CXX_STANDARD 14)

add_executable(test-CUASArgs
        CUASArgsTest.cpp
        )
target_link_libraries(test-CUASArgs cuascore)
target_link_libraries(test-CUASArgs gtest_main)
set_property(TARGET test-CUASArgs PROPERTY CXX_STANDARD 14)

add_executable(test-timeparse
        timeparseTest.cpp
        )
target_link_libraries(test-timeparse cuascore)
target_link_libraries(test-timeparse gtest_main)
set_property(TARGET test-timeparse PROPERTY CXX_STANDARD 14)
