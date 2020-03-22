set(gpds_source
  gpdevolutionkernels.cc
  gpdbuilder.cc
  operatorerbl.cc
  )

add_library(gpds OBJECT ${gpds_source})
