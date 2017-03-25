# build tests
add_executable(subgrid_test subgrid_test.cc)
target_link_libraries(subgrid_test apfel)
add_test(subgrid_test subgrid_test)

# build tests
add_executable(grid_test grid_test.cc)
target_link_libraries(grid_test apfel)
add_test(grid_test grid_test)

# build tests
add_executable(interpolation_test interpolation_test.cc)
target_link_libraries(interpolation_test apfel)
add_test(interpolation_test interpolation_test)

# build tests
add_executable(integrator_test integrator_test.cc)
target_link_libraries(integrator_test apfel)
add_test(integrator_test integrator_test)

# build tests
add_executable(operator_test operator_test.cc)
target_link_libraries(operator_test apfel)
add_test(operator_test operator_test)

# build tests
add_executable(distribution_test distribution_test.cc)
target_link_libraries(distribution_test apfel)
add_test(distribution_test distribution_test)

# build tests
add_executable(coupling_test coupling_test.cc)
target_link_libraries(coupling_test apfel)
add_test(coupling_test coupling_test)

# build tests
add_executable(structurefunction_test structurefunction_test.cc)
target_link_libraries(structurefunction_test apfel)
add_test(structurefunction_test structurefunction_test)

# build tests
add_executable(qgrid_test qgrid_test.cc)
target_link_libraries(qgrid_test apfel)
add_test(qgrid_test qgrid_test)

# build tests
add_executable(set_test set_test.cc)
target_link_libraries(set_test apfel)
add_test(set_test set_test)

# build tests
add_executable(cache_test cache_test.cc)
target_link_libraries(cache_test apfel)
add_test(cache_test cache_test)

# build tests
add_executable(qgriddist_test qgriddist_test.cc)
target_link_libraries(qgriddist_test apfel)
add_test(qgriddist_test qgriddist_test)

# build tests
add_executable(dglap_test dglap_test.cc)
target_link_libraries(dglap_test apfel)
add_test(dglap_test dglap_test)
