set(evolution_module_source
  alphaqcd.cc
  gridalphaqcd.cc
  dglap.cc
  )

add_library(evolution_module OBJECT ${evolution_module_source})
