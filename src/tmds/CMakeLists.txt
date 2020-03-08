set(tmds_source
  tmdbuilder.cc
  gammak.cc
  gammaf.cc
  kcs.cc
  matchingfunctionspdf.cc
  matchingfunctionsff.cc
  ogataquadrature.cc
  doubleexponentialquadrature.cc
  twobodyphasespace.cc
  )

add_library(tmds OBJECT ${tmds_source})
