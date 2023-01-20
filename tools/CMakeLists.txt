add_executable(cuas.exe main.cpp)
target_link_libraries(cuas.exe cuascore)
target_link_libraries(cuas.exe petscwrapper)

install(TARGETS cuas.exe DESTINATION bin)
