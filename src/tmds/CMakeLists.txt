set(tmds_source
  tmdbuilder.cc
  gammacusp.cc
  gammav.cc
  gammacs.cc
  matchingfunctionspdf.cc
  matchingfunctionsff.cc
  ogataquadrature.cc
  tmdcrosssections.cc
  hardfactors.cc
  )

add_library(tmds OBJECT ${tmds_source})
