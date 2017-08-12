set(dis_source
  structurefunctionbuilder.cc
  zeromasscoefficientfunctions.cc
  massivecoefficientfunctions.cc
  hqcoef.f
  massivezerocoefficientfunctions.cc
  )

add_library(dis OBJECT ${dis_source})
