set(apfel_source
  subgrid.cc
  grid.cc
  utils.cc
  splittings.cc
  opset.cc
  evolinit.cc
  evolsetup.cc
  distset.cc
)

#add_library(apfel-shared SHARED ${apfel_source})
add_library(apfel-static STATIC ${apfel_source})
#set_target_properties(apfel-shared PROPERTIES OUTPUT_NAME apfel CLEAN_DIRECT_OUTPUT 1)
set_target_properties(apfel-static PROPERTIES OUTPUT_NAME apfel CLEAN_DIRECT_OUTPUT 1)
