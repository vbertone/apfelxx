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

# build tests
add_executable(qgrid_test qgrid_test.cc)
target_link_libraries(qgrid_test apfel)
add_test(performance_qgrid_test qgrid_test)
add_test(memcheck_qgrid_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/qgrid_test)

# build tests
add_executable(set_test set_test.cc)
target_link_libraries(set_test apfel)
add_test(performance_set_test set_test)
add_test(memcheck_set_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/set_test)

# build tests
add_executable(cache_test cache_test.cc)
target_link_libraries(cache_test apfel)
add_test(performance_cache_test cache_test)
add_test(memcheck_cache_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/cache_test)

# build tests
add_executable(qgriddist_test qgriddist_test.cc)
target_link_libraries(qgriddist_test apfel)
add_test(performance_qgriddist_test qgriddist_test)
add_test(memcheck_qgriddist_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/qgriddist_test)

# build tests
add_executable(dglapqcd_test dglapqcd_test.cc)
target_link_libraries(dglapqcd_test apfel)
add_test(performance_dglapqcd_test dglapqcd_test)
add_test(memcheck_dglapqcd_test ${CMAKE_CURRENT_BINARY_DIR}/memtest.py ${CMAKE_CURRENT_BINARY_DIR}/dglapqcd_test)
