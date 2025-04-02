add_library(netcdfwrapper src/NetCDFFile.cpp)

target_include_directories(netcdfwrapper
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

add_netcdf(netcdfwrapper)

target_link_libraries(netcdfwrapper petscwrapper)

add_spdlog(netcdfwrapper)

install(TARGETS netcdfwrapper DESTINATION lib)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.h")
