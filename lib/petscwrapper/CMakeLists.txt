add_library(petscwrapper src/PETScGrid.cpp)

target_include_directories(petscwrapper
                           PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(petscwrapper ${PETSC_LIBRARIES})
target_include_directories(petscwrapper PUBLIC ${PETSC_INCLUDES})

add_spdlog(petscwrapper)

install(TARGETS petscwrapper DESTINATION lib)
install(
  DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
  DESTINATION include
  FILES_MATCHING
  PATTERN "*.h")
