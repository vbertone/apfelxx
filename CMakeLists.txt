cmake_minimum_required (VERSION 2.6)
# Set a default build type if none was specified
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to 'Release' as none was specified.")
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release"
    "MinSizeRel" "RelWithDebInfo")
endif()

if (NOT DEFINED CMAKE_MACOSX_RPATH)
   set(CMAKE_MACOSX_RPATH 0)
endif()

# the project name
project(apfel CXX Fortran)

# activating some global properties for the project
set(apfel_VERSION 4.1.0)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -fPIC")
#set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -Wall -Wextra -fPIC")
set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS}  -fPIC")
set(CMAKE_SHARED_LINKER_FLAGS -w)

# export version to file
configure_file(
  "${PROJECT_SOURCE_DIR}/inc/apfel/version.h.in"
  "${PROJECT_BINARY_DIR}/inc/apfel/version.h"
)

# Configuration script
set(prefix ${CMAKE_INSTALL_PREFIX})
set(exec_prefix "${prefix}")
set(includedir "${prefix}/include")
set(libdir "${prefix}/lib")
configure_file(
  "${PROJECT_SOURCE_DIR}/bin/apfelxx-config.in"
  "${PROJECT_SOURCE_DIR}/bin/apfelxx-config"
)

include_directories(${PROJECT_SOURCE_DIR}/inc)

# build target library
add_subdirectory(src)

option(BUILD_TESTS "Build tests code" ON)
if(BUILD_TESTS)
  enable_testing()
  add_subdirectory(tests)
endif(BUILD_TESTS)

install(FILES ${PROJECT_SOURCE_DIR}/bin/apfelxx-config
	   DESTINATION bin
	   PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE)
