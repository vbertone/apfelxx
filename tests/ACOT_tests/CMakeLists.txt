# build ACOT tests
file(GLOB testcodes *.cc )
foreach(testsource ${testcodes})
  GET_FILENAME_COMPONENT(filename ${testsource} NAME_WE)
  add_executable(${filename} ${filename}.cc)
  target_link_libraries(${filename} apfelxx)
  # add_test(${filename} ${filename})
endforeach()