add_library(cuascore
)

target_include_directories(cuascore PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(cuascore ${PETSC_LIBRARIES})
target_include_directories(cuascore PRIVATE ${PETSC_INCLUDES})

