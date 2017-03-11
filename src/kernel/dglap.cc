//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/dglap.h"
#include "apfel/ode.h"

#include <iostream>
#include <cmath>

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  Dglap::Dglap(function<Set<Operator>(int,double)>      const& SplittingFunctions,
	       function<Set<Operator>(bool,int,double)> const& MatchingConditions,
	       Set<Distribution>                        const& ObjRef,
	       double                                   const& MuDistRef,
	       vector<double>                           const& Masses,
	       vector<double>                           const& Thresholds,
	       int                                      const& nstep):
    MatchedEvolution(ObjRef, MuDistRef, Masses, Thresholds, nstep),
    _SplittingFunctions(SplittingFunctions),
    _MatchingConditions(MatchingConditions),
    _nstep(nstep)
  {
  }

  //_________________________________________________________________________________
  Dglap::Dglap(function<Set<Operator>(int,double)>      const& SplittingFunctions,
	       function<Set<Operator>(bool,int,double)> const& MatchingConditions,
	       Set<Distribution>                        const& ObjRef,
	       double                                   const& MuDistRef,
	       vector<double>                           const& Masses,
	       int                                      const& nstep):
    Dglap(SplittingFunctions, MatchingConditions, ObjRef, MuDistRef, Masses, Masses, nstep)
  {
  }

  //_________________________________________________________________________________
  Set<Distribution> Dglap::MatchObject(bool const& Up, int const& nf, Set<Distribution> const& f) const
  {
    // Get matching conditions
    auto MC = _MatchingConditions(Up, nf, _LogTh2M2[nf]);

    // Create the object 'g' with the same convolution map of the matching conditions
    // but containing the same objects of the input set of functions 'f'.
    Set<Distribution> g{MC.GetMap(),f.GetObjects()};

    // Convolute 'MC' and 'g'
    auto MO = MC * g;

    // Return the convoluted object with the map on the next evolution step
    return Set<Distribution>{_SplittingFunctions(nf+1, 0).GetMap(), MO.GetObjects()};
  }

  //_________________________________________________________________________________
  Set<Distribution> Dglap::Derivative(int const& nf, double const& t, Set<Distribution> const& f) const
  {
    return _SplittingFunctions(nf, exp(t/2)) * f;
  }

}

