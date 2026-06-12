# build tests
file(GLOB testcodes *.cc )
foreach(testsource ${testcodes})
  GET_FILENAME_COMPONENT(filename ${testsource} NAME_WE)
  add_executable(${filename} ${filename}.cc)
  target_link_libraries(${filename} apfelxx)
  # Each test runs through regression.py, which executes the binary once,
  # checks its exit code, and compares the filtered output against the
  # stored benchmark in tests/benchmarks/
  add_test(NAME ${filename}
    COMMAND python3 ${CMAKE_CURRENT_SOURCE_DIR}/regression.py
            ${CMAKE_CURRENT_BINARY_DIR}/${filename}
            ${CMAKE_CURRENT_SOURCE_DIR}/benchmarks/${filename}.out
  )
  set_tests_properties(${filename} PROPERTIES
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    LABELS "regression"
  )
endforeach()
