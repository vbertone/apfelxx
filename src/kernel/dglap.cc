//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/dglap.h"
#include "apfel/ode.h"
#include "apfel/tools.h"

#include <iostream>
#include <cmath>

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  Dglap::Dglap(function<Set<Operator>(int const&,double const&)>      const& SplittingFunctions,
	       function<Set<Operator>(bool,int const&,double const&)> const& MatchingConditions,
	       Set<Distribution>                                      const& ObjRef,
	       double                                                 const& MuDistRef,
	       vector<double>                                         const& Masses,
	       vector<double>                                         const& Thresholds,
	       int                                                    const& nsteps):
    MatchedEvolution(ObjRef, MuDistRef, Masses, Thresholds, nsteps),
    _SplittingFunctions(SplittingFunctions),
    _MatchingConditions(MatchingConditions)
  {
  }

  //_________________________________________________________________________________
  Dglap::Dglap(function<Set<Operator>(int const&,double const&)>      const& SplittingFunctions,
	       function<Set<Operator>(bool,int const&,double const&)> const& MatchingConditions,
	       Set<Distribution>                                      const& ObjRef,
	       double                                                 const& MuDistRef,
	       vector<double>                                         const& Masses,
	       int                                                    const& nsteps):
    Dglap(SplittingFunctions, MatchingConditions, ObjRef, MuDistRef, Masses, Masses, nsteps)
  {
  }

  //_________________________________________________________________________________
  Set<Distribution> Dglap::MatchObject(bool const& Up, int const& nf, Set<Distribution> const& f) const
  {
    // Get matching conditions.
    auto MC = _MatchingConditions(Up, nf, _LogTh2M2[nf]);

    // Create the object 'g' with the same convolution map of the
    // matching conditions but containing the same objects of the
    // input set of functions 'f'.
    Set<Distribution> g{MC.GetMap(),f.GetObjects()};

    // Convolute 'MC' and 'g'.
    auto MO = MC * g;

    // Return the convoluted object with the map on the next evolution
    // step.
    return Set<Distribution>{_SplittingFunctions((Up ? nf+1 : nf-1), 0).GetMap(), MO.GetObjects()};
  }

  //_________________________________________________________________________________
  Set<Distribution> Dglap::Derivative(int const& nf, double const& t, Set<Distribution> const& f) const
  {
    return _SplittingFunctions(nf, exp(t/2)) * f;
  }


  //_________________________________________________________________________________
  void Dglap::SetInitialDistributions(function<double(int const&, double const&)> const& InDistFunc)
  {
    // Compute number of active flavours the the PDF initial scale.
    int nf0 = NF(_MuRef, _Thresholds);

    // Allocate initial scale distributions.
    unordered_map<int,Distribution> DistMap;
    for (int i = 0; i <= 12; i++)
      DistMap.insert({i,Distribution{_ObjRef.at(0).GetGrid(), InDistFunc, i}});

    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution>{_SplittingFunctions(nf0, 0).GetMap(), DistMap});
  }

  //_________________________________________________________________________________
  void Dglap::SetInitialDistributions(function<unordered_map<int,double>(double const&)> const& InDistFunc)
  {
    // Compute number of active flavours the the PDF initial scale.
    int nf0 = NF(_MuRef, _Thresholds);

    // Allocate initial scale distributions.
    unordered_map<int,Distribution> DistMap = DistributionMap(_ObjRef.at(0).GetGrid(), InDistFunc);

    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution>{_SplittingFunctions(nf0, 0).GetMap(), DistMap});
  }

}

