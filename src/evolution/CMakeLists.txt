set(evolution_source
  alphaqcd.cc
  betaqcd.cc
  alphaqed.cc
  betaqed.cc
  dglapbuilder.cc
  splittingfunctions.cc
  matchingconditions.cc
  splittingfunctions_tl.cc
  matchingconditions_tl.cc
  lhtoypdfs.cc
  )

add_library(evolution OBJECT ${evolution_source})
