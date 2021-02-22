add_library(petscwrapper
  src/PETScVec.cpp
  src/PETScMat.cpp
  src/PETScGrid.cpp
)

target_include_directories(petscwrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(petscwrapper ${PETSC_LIBRARIES})
target_include_directories(petscwrapper PUBLIC ${PETSC_INCLUDES})
