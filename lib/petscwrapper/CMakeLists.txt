add_library(petscwrapper
  src/PetscVec.cpp
  src/PetscMat.cpp
  src/PetscGrid.cpp
  src/PetscSolver.cpp
)

target_include_directories(petscwrapper PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(petscwrapper ${PETSC_LIBRARIES})
target_include_directories(petscwrapper PRIVATE ${PETSC_INCLUDES})
