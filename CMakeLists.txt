# BASIC DEFINITIONS ========================================================================

# define minimum version of cmake
cmake_minimum_required (VERSION 3.5)

# define project name and its language
project(apfel CXX Fortran)

# define c++ standard and issue all the warning demanded by this standard
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fPIC")

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# VERSION ==============================================================================

set(apfelxx_VERSION 4.7.0)

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/inc/apfel/version.h.in"
  "${PROJECT_SOURCE_DIR}/inc/apfel/version.h"
)

# CONFIG SCRIPT =========================================================================

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file("${PROJECT_SOURCE_DIR}/bin/apfelxx-config.in" "${PROJECT_SOURCE_DIR}/bin/apfelxx-config")

# FINALIZE ==================================================================================

# include directory
include_directories(${PROJECT_SOURCE_DIR}/inc)

# generate list of source files
file(GLOB_RECURSE  source_files src/*)

# build test codes
enable_testing()
add_subdirectory(tests)

# define target library
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_library(apfelxx SHARED ${source_files})

# define libraries to be linked
target_link_libraries(apfelxx)

# install
install(TARGETS apfelxx
           RUNTIME DESTINATION bin
           LIBRARY DESTINATION lib)
install(DIRECTORY inc/apfel DESTINATION include)
install(FILES ${PROJECT_SOURCE_DIR}/bin/apfelxx-config DESTINATION bin
PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
