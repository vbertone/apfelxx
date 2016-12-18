# build tests
add_executable(test1 test1.cc)
target_link_libraries(test1 apfel)
add_test(performance_test1 test1)
add_test(memcheck_test1 ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/test1)
