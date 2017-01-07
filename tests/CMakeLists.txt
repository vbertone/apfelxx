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
add_test(performance_interpolation_test interpolation_test)
add_test(memcheck_interpolation_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/interpolation_test)

# build tests
add_executable(integrator_test integrator_test.cc)
target_link_libraries(integrator_test apfel)
add_test(performance_integrator_test integrator_test)
add_test(memcheck_integrator_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/integrator_test)

# build tests
add_executable(operator_test operator_test.cc)
target_link_libraries(operator_test apfel)
add_test(performance_operator_test operator_test)
add_test(memcheck_operator_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/operator_test)

# build tests
add_executable(distribution_test distribution_test.cc)
target_link_libraries(distribution_test apfel)
add_test(performance_distribution_test distribution_test)
add_test(memcheck_distribution_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/distribution_test)

# build tests
add_executable(coupling_test coupling_test.cc)
target_link_libraries(coupling_test apfel)
add_test(performance_coupling_test coupling_test)
add_test(memcheck_coupling_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/coupling_test)

# build tests
add_executable(structurefunction_test structurefunction_test.cc)
target_link_libraries(structurefunction_test apfel)
add_test(performance_structurefunction_test structurefunction_test)
add_test(memcheck_structurefunction_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/structurefunction_test)
