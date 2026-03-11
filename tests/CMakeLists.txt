# build tests
file(GLOB testcodes *.cc )
foreach(testsource ${testcodes})
  GET_FILENAME_COMPONENT(filename ${testsource} NAME_WE)
  add_executable(${filename} ${filename}.cc)
  target_link_libraries(${filename} apfelxx)
  add_test(${filename} ${CMAKE_CURRENT_BINARY_DIR}/${filename})
  set_tests_properties(${filename} PROPERTIES WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})
endforeach()
