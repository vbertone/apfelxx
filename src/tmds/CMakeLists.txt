set(tmds_source
  tmdbuilder.cc
  matchingfunctions.cc
  ogataquadrature.cc
  tmdcrosssections.cc
  )

add_library(tmds OBJECT ${tmds_source})
