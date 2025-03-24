add_executable(test-yamlcpp general.cpp)
target_link_libraries(test-yamlcpp gtest)
add_yamlcpp(test-yamlcpp)

install(TARGETS test-yamlcpp DESTINATION test/thirdparty)
install(FILES testconfig.yaml DESTINATION test/thirdparty)
