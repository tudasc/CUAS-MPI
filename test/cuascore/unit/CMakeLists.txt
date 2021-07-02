add_executable(test-systemmatrix
        systemmatrixTest.cpp
)
target_link_libraries(test-systemmatrix cuascore)
target_link_libraries(test-systemmatrix gtest)
set_property(TARGET test-systemmatrix PROPERTY CXX_STANDARD 14)

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

add_executable(test-forcing
        forcingTest.cpp
        )
target_link_libraries(test-forcing cuascore)
target_link_libraries(test-forcing gtest)
set_property(TARGET test-forcing PROPERTY CXX_STANDARD 14)

add_executable(test-netcdf
        netcdfTest.cpp
        )
target_link_libraries(test-netcdf cuascore)
target_link_libraries(test-netcdf gtest)
set_property(TARGET test-netcdf PROPERTY CXX_STANDARD 14)

install(TARGETS test-systemmatrix DESTINATION test/cuascore/unit)
install(TARGETS test-CUASKernels DESTINATION test/cuascore/unit)
install(TARGETS test-specialgradient DESTINATION test/cuascore/unit)
install(TARGETS test-fillNoData DESTINATION test/cuascore/unit)
install(TARGETS test-CUASArgs DESTINATION test/cuascore/unit)
install(TARGETS test-timeparse DESTINATION test/cuascore/unit)
install(TARGETS test-forcing DESTINATION test/cuascore/unit)
install(TARGETS test-netcdf DESTINATION test/cuascore/unit)

