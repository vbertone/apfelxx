set(apfel_source
  subgrid.cc
  grid.cc
  interpolator.cc
  lagrangeinterpolator.cc
  integrator.cc
  operator.cc
  expression.cc
  #utils.cc
  #splittings.cc
  #opset.cc
  #evolinit.cc
  #evolsetup.cc
  #distset.cc
)

option(SHARED "Build shared-libray instead of static-libray" OFF)
if(SHARED)
  add_library(apfel SHARED ${apfel_source})
else(SHARED)
  add_library(apfel STATIC ${apfel_source})
endif(SHARED)
