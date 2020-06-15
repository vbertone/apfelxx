# build tests
add_executable(subgrid_test subgrid_test.cc)
target_link_libraries(subgrid_test apfelxx)
add_test(subgrid_test subgrid_test)

# build tests
add_executable(grid_test grid_test.cc)
target_link_libraries(grid_test apfelxx)
add_test(grid_test grid_test)

# build tests
add_executable(interpolation_test interpolation_test.cc)
target_link_libraries(interpolation_test apfelxx)
add_test(interpolation_test interpolation_test)

# build tests
add_executable(integrator_test integrator_test.cc)
target_link_libraries(integrator_test apfelxx)
add_test(integrator_test integrator_test)

# build tests
add_executable(operator_test operator_test.cc)
target_link_libraries(operator_test apfelxx)
add_test(operator_test operator_test)

# build tests
add_executable(distribution_test distribution_test.cc)
target_link_libraries(distribution_test apfelxx)
add_test(distribution_test distribution_test)

# build tests
add_executable(coupling_test coupling_test.cc)
target_link_libraries(coupling_test apfelxx)
add_test(coupling_test coupling_test)

# build tests
add_executable(structurefunction_test structurefunction_test.cc)
target_link_libraries(structurefunction_test apfelxx)
add_test(structurefunction_test structurefunction_test)

# build tests
add_executable(qgrid_test qgrid_test.cc)
target_link_libraries(qgrid_test apfelxx)
add_test(qgrid_test qgrid_test)

# build tests
add_executable(set_test set_test.cc)
target_link_libraries(set_test apfelxx)
add_test(set_test set_test)

# build tests
add_executable(qgriddist_test qgriddist_test.cc)
target_link_libraries(qgriddist_test apfelxx)
add_test(qgriddist_test qgriddist_test)

# build tests
add_executable(dglap_test dglap_test.cc)
target_link_libraries(dglap_test apfelxx)
add_test(dglap_test dglap_test)

# build tests
add_executable(ogata_test ogata_test.cc)
target_link_libraries(ogata_test apfelxx)
add_test(ogata_test ogata_test)

# build tests
add_executable(gpd_test gpd_test.cc)
target_link_libraries(gpd_test apfelxx)
add_test(gpd_test gpd_test)

# build tests
add_executable(polylog_test polylog_test.cc)
target_link_libraries(polylog_test apfelxx)
add_test(polylog_test polylog_test)

# build tests
add_executable(expansion_test expansion_test.cc)
target_link_libraries(expansion_test apfelxx)
add_test(expansion_test expansion_test)

# build tests
add_executable(siaxsec_test siaxsec_test.cc)
target_link_libraries(siaxsec_test apfelxx)
add_test(siaxsec_test siaxsec_test)
