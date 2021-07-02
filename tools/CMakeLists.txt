add_executable(cuas.exe
	main.cpp
)
target_link_libraries(cuas.exe cuascore)
target_link_libraries(cuas.exe petscwrapper) 
set_property(TARGET cuas.exe PROPERTY CXX_STANDARD 14)

install(TARGETS cuas.exe DESTINATION bin)
