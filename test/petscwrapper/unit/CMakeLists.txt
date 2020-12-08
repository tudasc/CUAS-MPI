add_executable(vec-test
  VecTest.cpp
)
target_link_libraries(vec-test cuascore)

add_executable(mat-test
  MatTest.cpp
)
target_link_libraries(mat-test cuascore)

add_executable(grid-test
  GridTest.cpp
)
target_link_libraries(grid-test cuascore)
