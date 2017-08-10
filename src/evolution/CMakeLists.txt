set(evolution_source
  alphaqcd.cc
  dglapbuilder.cc
  splittingfunctions.cc
  matchingconditions.cc
  )

add_library(evolution OBJECT ${evolution_source})
