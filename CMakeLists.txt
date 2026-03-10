# ============================
# BASIC DEFINITIONS
# ============================

# Define minimum version of cmake
cmake_minimum_required(VERSION 3.10)

# Define project name and its language
project(apfel CXX Fortran)

# Define code version
set(apfelxx_VERSION 4.10.0)

# Define c++ standard and issue all the warning demanded by this standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC -O3 -Wno-psabi")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fPIC -O3")

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif ()

if (CMAKE_VERSION VERSION_GREATER_EQUAL "3.27.0")
   cmake_policy(SET CMP0148 OLD)
endif ()

# ============================
# MAIN
# ============================

# Define include directory
include_directories(${PROJECT_SOURCE_DIR}/inc)

# Generate list of source files
file(GLOB_RECURSE  source_files src/*)

# Define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
add_library(apfelxx SHARED ${source_files})

target_include_directories(apfelxx
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<INSTALL_INTERFACE:include>
)

# Define libraries to be linked
target_link_libraries(apfelxx PUBLIC)

# LHAPDF
find_program(LHAPDF_CONFIG lhapdf-config)

if (LHAPDF_CONFIG)
   message ("-- LHAPDF found. The 'LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const&)' function will be built.")
   set(WITH_LHAPDF 1)
   execute_process(COMMAND ${LHAPDF_CONFIG} --cflags OUTPUT_VARIABLE LHAPDF_CXX_FLAGS OUTPUT_STRIP_TRAILING_WHITESPACE)
   execute_process(COMMAND ${LHAPDF_CONFIG} --libs OUTPUT_VARIABLE LHAPDF_LIBRARIES OUTPUT_STRIP_TRAILING_WHITESPACE)
   set(LHAPDF_CXX_FLAGS ${LHAPDF_CXX_FLAGS} CACHE STRING INTERNAL)
   set(LHAPDF_LIBRARIES ${LHAPDF_LIBRARIES} CACHE STRING INTERNAL)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${LHAPDF_CXX_FLAGS}")
   target_link_libraries(apfelxx PUBLIC ${LHAPDF_LIBRARIES})
else (LHAPDF_CONFIG)
   message ("-- LHAPDF not found. The 'LHAPDF::PDF* mkPDF(apfel::InitialiseEvolution const&)' function will not be built.")
endif (LHAPDF_CONFIG)


find_package(yaml-cpp QUIET)
if (yaml-cpp_FOUND)
   message ("   + yaml-cpp found. DoubleExpression objects can be exported to file.")
   set(WITH_YAML_CPP 1)
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${YAML_CPP_CFLAGS}")
   target_link_libraries(apfelxx PUBLIC yaml-cpp::yaml-cpp)
else (yaml-cpp_FOUND)
   message ("   + yaml-cpp not found. DoubleExpression objects cannot be exported to file.")
endif (yaml-cpp_FOUND)


# pybind11
find_package(pybind11 QUIET)

if (pybind11_FOUND)
   message ("-- pybind11 found. The python wrapper will be built.")
   pybind11_add_module(apfelpy pywrap/pywrap.cc)
   target_link_libraries(apfelpy PRIVATE apfelxx)
   execute_process(COMMAND "${PYTHON_EXECUTABLE}" -c "import sysconfig; print('${CMAKE_INSTALL_PREFIX}' + sysconfig.get_path('purelib').replace(sysconfig.get_path('data'), ''), end = '')" OUTPUT_VARIABLE PYTHON_INSTALLATION_PATH)
   install(TARGETS apfelpy LIBRARY DESTINATION "${PYTHON_INSTALLATION_PATH}")
else (pybind11_FOUND)
   message ("-- pybind11 not found. The python wrapper will not be built.")
endif (pybind11_FOUND)

# Enable testing and build test codes
enable_testing()
add_subdirectory(tests)

# ============================
# CONFIGURATION
# ============================

# Configure version to file
configure_file("${PROJECT_SOURCE_DIR}/inc/apfel/version.h.in" "${PROJECT_SOURCE_DIR}/inc/apfel/version.h")

# Configure configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
set(exec_prefix "${prefix}")
set(cflags "-I${prefix}/include")
set(ldflags "-L${prefix}/lib -lapfelxx")
if (LHAPDF_CONFIG)
  set(cflags "${cflags} ${LHAPDF_CXX_FLAGS}")
  set(ldflags "${ldflags} ${LHAPDF_LIBRARIES}")
endif (LHAPDF_CONFIG)
if (YAML_CPP_FOUND)
  set(cflags "${cflags} ${YAML_CPP_CFLAGS}")
  set(ldflags "${ldflags} -L${YAML_CPP_LIBRARY_DIRS} -l${YAML_CPP_LIBRARIES}")
endif (YAML_CPP_FOUND)

# Configuration script
configure_file("${PROJECT_SOURCE_DIR}/bin/apfelxx-config.in" "${PROJECT_BINARY_DIR}/bin/apfelxx-config")

# Configure config.h header file
configure_file("${PROJECT_SOURCE_DIR}/inc/apfel/config.h.in" "${PROJECT_SOURCE_DIR}/inc/apfel/config.h")

# pkg-config
configure_file("${PROJECT_SOURCE_DIR}/bin/apfelxx.pc.in" "${PROJECT_BINARY_DIR}/bin/apfelxx.pc")

# Config.cmake file
include(CMakePackageConfigHelpers)
configure_package_config_file(${CMAKE_CURRENT_SOURCE_DIR}/cmake/apfelxxConfig.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfig.cmake INSTALL_DESTINATION lib/cmake/apfelxx)
write_basic_package_version_file(${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfigVersion.cmake VERSION ${apfelxx_VERSION} COMPATIBILITY SameMajorVersion)

# ============================
# INSTALLATION
# ============================

# Installation of code
install(TARGETS apfelxx EXPORT apfelxxTargets RUNTIME DESTINATION bin LIBRARY DESTINATION lib ARCHIVE DESTINATION lib)
install(DIRECTORY inc/apfel DESTINATION include)
install(FILES ${PROJECT_BINARY_DIR}/bin/apfelxx-config DESTINATION bin PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
install(FILES ${PROJECT_BINARY_DIR}/bin/apfelxx.pc DESTINATION lib/pkgconfig)
install(EXPORT apfelxxTargets FILE apfelxxTargets.cmake NAMESPACE apfelxx:: DESTINATION lib/cmake/apfelxx)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfig.cmake ${CMAKE_CURRENT_BINARY_DIR}/apfelxxConfigVersion.cmake DESTINATION lib/cmake/apfelxx)
