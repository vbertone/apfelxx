add_subdirectory(kernel)
add_subdirectory(evolution)
add_subdirectory(dis)

option(SHARED "Build shared-libray instead of static-libray" OFF)
if(SHARED)
  add_library(apfelxx SHARED $<TARGET_OBJECTS:kernel>
  $<TARGET_OBJECTS:evolution> $<TARGET_OBJECTS:dis>)
else(SHARED)
  add_library(apfelxx STATIC $<TARGET_OBJECTS:kernel>
  $<TARGET_OBJECTS:evolution> $<TARGET_OBJECTS:dis>)
endif(SHARED)

install(DIRECTORY ${PROJECT_SOURCE_DIR}/inc/apfel DESTINATION include)
install(TARGETS apfelxx DESTINATION lib)
