set(dis_source
  structurefunctionbuilder.cc
  zeromasscoefficientfunctions.cc
  massivecoefficientfunctions.cc
  hqcoef.f
  )

add_library(dis OBJECT ${dis_source})
