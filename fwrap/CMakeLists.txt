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

# Propagate the .mod search path to anything that links this target, both
# from within the build tree (other targets in this directory already do
# this explicitly below) and after installation.
target_include_directories(apfelxxfortran
    INTERFACE
        $<BUILD_INTERFACE:${CMAKE_Fortran_MODULE_DIRECTORY}>
        $<INSTALL_INTERFACE:include/apfel/fortran>
)

# .mod files are compiler-specific binary interface files, unlike the
# portable C++ headers under include/apfel, so they get their own
# subdirectory rather than being installed alongside them.
install(TARGETS apfelxxfortran EXPORT apfelxxTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY ${CMAKE_Fortran_MODULE_DIRECTORY}/ DESTINATION include/apfel/fortran)

add_executable(checkapfel_test checkapfel_test.f90)
target_link_libraries(checkapfel_test PRIVATE apfelxxfortran)
target_include_directories(checkapfel_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(initialize_test initialize_test.f90)
target_link_libraries(initialize_test PRIVATE apfelxxfortran)
target_include_directories(initialize_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(evolve_test evolve_test.f90)
target_link_libraries(evolve_test PRIVATE apfelxxfortran)
target_include_directories(evolve_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(renfacratio_test renfacratio_test.f90)
target_link_libraries(renfacratio_test PRIVATE apfelxxfortran)
target_include_directories(renfacratio_test PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})

add_executable(tabulation tabulation.f90)
target_link_libraries(tabulation PRIVATE apfelxxfortran)
target_include_directories(tabulation PRIVATE ${CMAKE_Fortran_MODULE_DIRECTORY})
