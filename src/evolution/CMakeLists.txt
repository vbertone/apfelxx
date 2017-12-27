set(evolution_source
  alphaqcd.cc
  betaqcd.cc
  dglapbuilder.cc
  splittingfunctions.cc
  matchingconditions.cc
  splittingfunctions_tl.cc
  matchingconditions_tl.cc
  )

add_library(evolution OBJECT ${evolution_source})
