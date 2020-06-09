set(structurefunctions_source
  structurefunctionbuilder.cc
  zeromasscoefficientfunctions.cc
  massivecoefficientfunctions.cc
  hqcoef.f
  massivezerocoefficientfunctions.cc
  zeromasscoefficientfunctions_tl.cc
  disbasis.cc
  )

add_library(structurefunctions OBJECT ${structurefunctions_source})

set_property(TARGET structurefunctions PROPERTY POSITION_INDEPENDENT_CODE ON)
