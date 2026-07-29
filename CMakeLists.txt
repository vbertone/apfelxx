# ============================
# BASIC DEFINITIONS
# ============================

cmake_minimum_required(VERSION 3.10)
project(apfel CXX Fortran)

set(apfelxx_VERSION 4.10.0)

# Build type and C++ standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Enable position independent code globally
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Optional code-coverage instrumentation (gcov / llvm-cov). Default OFF, so a
# normal build is completely unaffected. When ON, instrument every C++ target
# (library, tests and the python wrapper) so that inline code from the headers
# in inc/ is measured as well. The -O0 -g flags come after the per-target -O3
# below and, since the last flag wins, correctly disable optimisation in
# coverage builds only. See tests/coverage.sh for the report-generation step.
option(APFEL_COVERAGE "Build with code-coverage instrumentation (gcov/llvm-cov)" OFF)
if(APFEL_COVERAGE)
    message("-- Code coverage enabled (instrumented build: --coverage -O0 -g).")
    add_compile_options(
        $<$<COMPILE_LANGUAGE:CXX>:--coverage>
        $<$<COMPILE_LANGUAGE:CXX>:-O0>
        $<$<COMPILE_LANGUAGE:CXX>:-g>
        $<$<COMPILE_LANGUAGE:CXX>:-fprofile-update=atomic>)
    add_link_options(--coverage)
endif()

# ============================
# MAIN
# ============================

# Generate list of source files
file(GLOB_RECURSE source_files src/*)

set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
add_library(apfelxx SHARED ${source_files})

# Target-specific include directories
target_include_directories(apfelxx
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

# Compile options per target instead of global flags
target_compile_options(apfelxx PRIVATE
    $<$<COMPILE_LANGUAGE:CXX>:-Wall -Wextra -O3 -Wno-psabi>
    $<$<COMPILE_LANGUAGE:Fortran>:-O3>
)

# Set MODULES path
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules ${CMAKE_MODULE_PATH})

# LHAPDF
find_package(LHAPDF QUIET)
if(LHAPDF_FOUND)
    message("-- LHAPDF found. The 'LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const&)' function will be built.")
    set(WITH_LHAPDF 1)
    target_link_libraries(apfelxx PUBLIC LHAPDF::LHAPDF)
    # No need to manually extract cflags/libs – CMake handles it.
else()
    message("-- LHAPDF not found. The 'LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const&)' function will not be built.")
endif()

# cereal: header-only and bundled under inc/apfel/cereal, so DoubleOperator
# binary (de)serialisation is always available with no external dependency.

# pybind11
set(PYBIND11_FINDPYTHON ON)
find_package(pybind11 QUIET)

if(pybind11_FOUND)
    message("-- pybind11 found. The python wrapper will be built.")
    pybind11_add_module(apfelpy pywrap/pywrap.cc)
    target_link_libraries(apfelpy PRIVATE apfelxx)

    # The installed module is placed in the Python site-packages directory,
    # away from libapfelxx (installed in ${CMAKE_INSTALL_PREFIX}/lib), whose
    # install name is @rpath/libapfelxx.dylib. Without an install RPATH the
    # build-tree RPATH is stripped on install and the loader cannot resolve
    # @rpath, so bake the library location into the installed module.
    set_target_properties(apfelpy PROPERTIES INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/lib")

    execute_process(
        COMMAND "${Python_EXECUTABLE}" -c
        "import sysconfig; print('${CMAKE_INSTALL_PREFIX}' + sysconfig.get_path('purelib').replace(sysconfig.get_path('data'), ''), end = '')"
        OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH
    )
    install(TARGETS apfelpy LIBRARY DESTINATION "${PYTHON_INSTALLATION_PATH}")
else()
    message("-- pybind11 not found. The python wrapper will not be built.")
endif()

# Testing
enable_testing()

# Fortran wrapper: mirrors the old Fortran-native APFEL interface
option(APFEL_FORTRAN_WRAPPER "Build the Fortran wrapper (mirrors the old Fortran-native APFEL interface)" ON)
if(APFEL_FORTRAN_WRAPPER)
    message("-- The Fortran wrapper will be built.")
    add_subdirectory(fwrap)
else()
    message("-- The Fortran wrapper will not be built.")
endif()

add_subdirectory(tests)

# ============================
# CONFIGURATION
# ============================

# Version header
configure_file(
    "${PROJECT_SOURCE_DIR}/inc/apfel/version.h.in"
    "${PROJECT_SOURCE_DIR}/inc/apfel/version.h"
)

# Helper script apfelxx-config
# Generate cflags/ldflags using target properties where possible.
# Here we keep the simple approach because it's used by the config script.
set(prefix ${CMAKE_INSTALL_PREFIX})
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
set(exec_prefix "${prefix}")
set(cflags "-I${prefix}/include")
set(ldflags "-L${prefix}/lib -lapfelxx")

if(LHAPDF_FOUND)
    # For LHAPDF, we can retrieve flags from the imported target.
    get_target_property(LHAPDF_INCLUDE_DIRS LHAPDF::LHAPDF INTERFACE_INCLUDE_DIRECTORIES)
    get_target_property(LHAPDF_LINK_LIBRARIES LHAPDF::LHAPDF INTERFACE_LINK_LIBRARIES)
    if(LHAPDF_INCLUDE_DIRS)
        string(REPLACE ";" " " LHAPDF_CXX_FLAGS "-I${LHAPDF_INCLUDE_DIRS}")
        set(cflags "${cflags} ${LHAPDF_CXX_FLAGS}")
    endif()
    if(LHAPDF_LINK_LIBRARIES)
        set(ldflags "${ldflags} ${LHAPDF_LINK_LIBRARIES}")
    endif()
endif()

# Fortran wrapper flags (empty if APFEL_FORTRAN_WRAPPER was OFF, so
# apfelxx-config's --fortran-cflags/--fortran-libs just print nothing
# rather than referencing a library that was never built).
set(fortran_cflags "")
set(fortran_ldflags "")
if(APFEL_FORTRAN_WRAPPER)
    set(fortran_cflags "-I${prefix}/include/apfel/fortran")
    set(fortran_ldflags "-L${prefix}/lib -lapfelxxfortran")
endif()

configure_file(
    "${PROJECT_SOURCE_DIR}/bin/apfelxx-config.in"
    "${PROJECT_BINARY_DIR}/bin/apfelxx-config"
)

# config.h header
configure_file(
    "${PROJECT_SOURCE_DIR}/inc/apfel/config.h.in"
    "${PROJECT_SOURCE_DIR}/inc/apfel/config.h"
)

# pkg-config
configure_file(
    "${PROJECT_SOURCE_DIR}/bin/apfelxx.pc.in"
    "${PROJECT_BINARY_DIR}/bin/apfelxx.pc"
)

# pkg-config for the Fortran wrapper: a separate .pc file (rather than
# folding it into apfelxx.pc), since it's an optional component with its
# own library and include directory; "Requires: apfelxx" pulls in the
# C++ library's own flags transitively.
if(APFEL_FORTRAN_WRAPPER)
    configure_file(
        "${PROJECT_SOURCE_DIR}/bin/apfelxxfortran.pc.in"
        "${PROJECT_BINARY_DIR}/bin/apfelxxfortran.pc"
    )
endif()

# CMake config files
include(CMakePackageConfigHelpers)
# INCLUDE_INSTALL_DIR feeds apfelxxConfig.cmake.in's @PACKAGE_INCLUDE_INSTALL_DIR@
# (via PATH_VARS below), which is what apfelxx_INCLUDE_DIRS resolves to.
set(INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_PREFIX}/include)
configure_package_config_file(
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/apfelxxConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfig.cmake
    INSTALL_DESTINATION lib/cmake/apfelxx
    PATH_VARS INCLUDE_INSTALL_DIR
)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfigVersion.cmake
    VERSION ${apfelxx_VERSION}
    COMPATIBILITY SameMajorVersion
)

# ============================
# INSTALLATION
# ============================

install(TARGETS apfelxx EXPORT apfelxxTargets
    RUNTIME DESTINATION bin
    LIBRARY DESTINATION lib
    ARCHIVE DESTINATION lib
)
install(DIRECTORY inc/apfel DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/bin/apfelxx-config
    DESTINATION bin
    PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)
install(FILES ${PROJECT_BINARY_DIR}/bin/apfelxx.pc
    DESTINATION lib/pkgconfig
)
if(APFEL_FORTRAN_WRAPPER)
    install(FILES ${PROJECT_BINARY_DIR}/bin/apfelxxfortran.pc
        DESTINATION lib/pkgconfig
    )
endif()
install(EXPORT apfelxxTargets
    FILE apfelxxTargets.cmake
    NAMESPACE apfelxx::
    DESTINATION lib/cmake/apfelxx
)
install(FILES
    ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfig.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfigVersion.cmake
    DESTINATION lib/cmake/apfelxx
)
install(FILES
    ${CMAKE_CURRENT_SOURCE_DIR}/cmake/Modules/FindLHAPDF.cmake
    DESTINATION lib/cmake/apfelxx
)
