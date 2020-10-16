add_library(cuascore
  src/PetscVec.cpp
  src/PetscMat.cpp
  src/PetscGrid.cpp
  src/specialgradient.cpp
  src/CUASModel.cpp
)

target_include_directories(cuascore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(cuascore ${PETSC_LIBRARIES})
target_include_directories(cuascore PRIVATE ${PETSC_INCLUDES})

