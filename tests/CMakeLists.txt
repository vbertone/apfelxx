# build tests
add_executable(subgrid_test subgrid_test.cc)
target_link_libraries(subgrid_test apfel)
add_test(performance_subgrid_test subgrid_test)
add_test(memcheck_subgrid_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/subgrid_test)
