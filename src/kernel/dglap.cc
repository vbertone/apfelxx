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
  template<class T>
  Dglap<T>::Dglap(function<Set<Operator>(int const&,double const&)> const& SplittingFunctions,
		  function<Set<Operator>(bool const&,int const&)>   const& MatchingConditions,
		  Set<T>                                            const& ObjRef,
		  double                                            const& MuDistRef,
		  vector<double>                                    const& Thresholds,
		  int                                               const& nsteps):
    MatchedEvolution<Set<T>>(ObjRef, MuDistRef, Thresholds, nsteps),
    _SplittingFunctions(SplittingFunctions),
    _MatchingConditions(MatchingConditions)
  {
  }

  //_________________________________________________________________________________
  template<class T>
  Set<T> Dglap<T>::MatchObject(bool const& Up, int const& nf, Set<T> const& f) const
  {
    // Get matching conditions.
    Set<Operator> MC = _MatchingConditions(Up, nf);

    // Create the object 'g' with the same convolution map of the
    // matching conditions but containing the same objects of the
    // input set of functions 'f'.
    Set<T> g{MC.GetMap(),f.GetObjects()};

    // Convolute 'MC' and 'g'.
    Set<T> MO = MC * g;
    MO.SetMap(_SplittingFunctions((Up ? nf+1 : nf-1), 0).GetMap());

    // Return the convoluted object with the map on the next evolution
    // step.
    return MO;
  }

  //_________________________________________________________________________________
  template<class T>
  Set<T> Dglap<T>::Derivative(int const& nf, double const& t, Set<T> const& f) const
  {
    return _SplittingFunctions(nf, exp(t/2)) * f;
  }

  // Fixed template types.
  template class Dglap<Distribution>;
  template class Dglap<Operator>;

  //_________________________________________________________________________________
  template<>
  void Dglap<Distribution>::SetInitialDistributions(function<double(int const&, double const&)> const& InDistFunc)
  {
    // Allocate initial scale distributions.
    map<int,Distribution> DistMap;
    for (int i = 0; i <= 12; i++)
      DistMap.insert({i,Distribution{_ObjRef.at(0).GetGrid(), InDistFunc, i}});

    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution>{_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistMap});
  }

  //_________________________________________________________________________________
  template<>
  void Dglap<Distribution>::SetInitialDistributions(function<map<int,double>(double const&)> const& InDistFunc)
  {
    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution>{_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistributionMap(_ObjRef.at(0).GetGrid(), InDistFunc)});
  }
}

