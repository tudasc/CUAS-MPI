add_executable(test-systemmatrix systemmatrixTest.cpp)
target_link_libraries(test-systemmatrix cuascore)
target_link_libraries(test-systemmatrix gtest)

add_executable(test-CUASKernels CUASKernelsTest.cpp)
target_link_libraries(test-CUASKernels cuascore)
target_link_libraries(test-CUASKernels gtest)

add_executable(test-fillNoData fillNoDataTest.cpp)
target_link_libraries(test-fillNoData cuascore)
target_link_libraries(test-fillNoData gtest)

add_executable(test-CUASArgs CUASArgsTest.cpp)
target_link_libraries(test-CUASArgs cuascore)
target_link_libraries(test-CUASArgs gtest)

add_executable(test-timeparse timeparseTest.cpp)
target_link_libraries(test-timeparse cuascore)
target_link_libraries(test-timeparse gtest)

add_executable(test-utilities utilitiesTest.cpp)
target_link_libraries(test-utilities cuascore)
target_link_libraries(test-utilities gtest)

add_executable(test-forcing forcingTest.cpp)
target_link_libraries(test-forcing cuascore)
target_link_libraries(test-forcing gtest)

add_executable(test-ModelReader ModelReaderTest.cpp)
target_link_libraries(test-ModelReader cuascore)
target_link_libraries(test-ModelReader gtest)

add_executable(test-SolutionHandler SolutionHandlerTest.cpp)
target_link_libraries(test-SolutionHandler cuascore)
target_link_libraries(test-SolutionHandler gtest)

add_executable(test-initialHead initialHeadTest.cpp)
target_link_libraries(test-initialHead cuascore)
target_link_libraries(test-initialHead gtest)

add_executable(test-timeForcingNetcdf TimeForcingFromNetcdfTest.cpp)
target_link_libraries(test-timeForcingNetcdf cuascore)
target_link_libraries(test-timeForcingNetcdf gtest)

add_executable(test-TimeIntegrator TimeIntegratorTest.cpp)
target_link_libraries(test-TimeIntegrator cuascore)
target_link_libraries(test-TimeIntegrator gtest)

install(TARGETS test-systemmatrix DESTINATION test/cuascore/unit)
install(TARGETS test-CUASKernels DESTINATION test/cuascore/unit)
install(TARGETS test-fillNoData DESTINATION test/cuascore/unit)
install(TARGETS test-CUASArgs DESTINATION test/cuascore/unit)
install(FILES testconfig.yaml DESTINATION test/cuascore/unit)
install(TARGETS test-timeparse DESTINATION test/cuascore/unit)
install(TARGETS test-utilities DESTINATION test/cuascore/unit)
install(TARGETS test-forcing DESTINATION test/cuascore/unit)
install(TARGETS test-ModelReader DESTINATION test/cuascore/unit)
install(TARGETS test-SolutionHandler DESTINATION test/cuascore/unit)
install(TARGETS test-initialHead DESTINATION test/cuascore/unit)
install(TARGETS test-timeForcingNetcdf DESTINATION test/cuascore/unit)
install(TARGETS test-TimeIntegrator DESTINATION test/cuascore/unit)
