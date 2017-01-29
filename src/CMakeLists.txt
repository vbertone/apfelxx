add_subdirectory(kernel)
add_subdirectory(evolution_module)

option(SHARED "Build shared-libray instead of static-libray" OFF)
if(SHARED)
  add_library(apfel SHARED $<TARGET_OBJECTS:kernel> $<TARGET_OBJECTS:evolution_module>)
else(SHARED)
  add_library(apfel STATIC $<TARGET_OBJECTS:kernel> $<TARGET_OBJECTS:evolution_module>)
endif(SHARED)
