set(tmds_source
  tmdbuilder.cc
  gammacusp.cc
  gammav.cc
  gammacs.cc
  logzeta.cc
  matchingfunctions.cc
  ogataquadrature.cc
  tmdcrosssections.cc
  )

add_library(tmds OBJECT ${tmds_source})
