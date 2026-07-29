# Fortran wrapper: mirrors the old Fortran-native APFEL interface on top
# of apfelxx. Both sides of the C++/Fortran boundary use standard
# extern "C" / bind(C) linkage with explicit, matching symbol names, so
# no compiler-specific name-mangling detection (e.g. CMake's
# FortranCInterface module) is needed.

set(CMAKE_Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/mod)

add_library(apfelxxfortran SHARED
    apfelfortran.cc
    apfel_fortran.f90
)
target_link_libraries(apfelxxfortran PRIVATE apfelxx)

add_executable(checkapfel_test checkapfel_test.f90)
target_link_libraries(checkapfel_test PRIVATE apfelxxfortran)
target_include_directories(checkapfel_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME checkapfel_test COMMAND checkapfel_test)

add_executable(initialize_test initialize_test.f90)
target_link_libraries(initialize_test PRIVATE apfelxxfortran)
target_include_directories(initialize_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME initialize_test COMMAND initialize_test)

add_executable(evolve_test evolve_test.f90)
target_link_libraries(evolve_test PRIVATE apfelxxfortran)
target_include_directories(evolve_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME evolve_test COMMAND evolve_test)

add_executable(renfacratio_test renfacratio_test.f90)
target_link_libraries(renfacratio_test PRIVATE apfelxxfortran)
target_include_directories(renfacratio_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
add_test(NAME renfacratio_test COMMAND renfacratio_test)

# Example program (mirrors the old examples/Tabulation.f), not a ctest:
# it reads Q0^2, Q^2 interactively from stdin, same as the original.
add_executable(tabulation tabulation.f90)
target_link_libraries(tabulation PRIVATE apfelxxfortran)
target_include_directories(tabulation PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
