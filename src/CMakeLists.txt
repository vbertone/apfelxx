add_subdirectory(kernel)
add_subdirectory(interface)
add_subdirectory(evolution)
add_subdirectory(dis)
add_subdirectory(tmds)

option(SHARED "Build shared-libray instead of static-libray" OFF)
if(SHARED)
  add_library(apfelxx SHARED $<TARGET_OBJECTS:kernel> $<TARGET_OBJECTS:interface>
  $<TARGET_OBJECTS:evolution> $<TARGET_OBJECTS:dis> $<TARGET_OBJECTS:tmds>)
else(SHARED)
  add_library(apfelxx STATIC $<TARGET_OBJECTS:kernel> $<TARGET_OBJECTS:interface>
  $<TARGET_OBJECTS:evolution> $<TARGET_OBJECTS:dis> $<TARGET_OBJECTS:tmds>)
endif(SHARED)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/inc/apfel DESTINATION include)
install(TARGETS apfelxx DESTINATION lib)
