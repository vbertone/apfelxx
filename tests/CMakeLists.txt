# build tests
add_executable(subgrid_test subgrid_test.cc)
target_link_libraries(subgrid_test apfel)
add_test(performance_subgrid_test subgrid_test)
add_test(memcheck_subgrid_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/subgrid_test)

# build tests
add_executable(grid_test grid_test.cc)
target_link_libraries(grid_test apfel)
add_test(performance_grid_test grid_test)
add_test(memcheck_grid_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/grid_test)

# build tests
add_executable(interpolation_test interpolation_test.cc)
target_link_libraries(interpolation_test apfel)
add_test(performance_interpolation_test grid_test)
add_test(memcheck_interpolation_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/interpolation_test)
