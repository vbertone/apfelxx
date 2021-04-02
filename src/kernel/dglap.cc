//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/dglap.h"
#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________________
  template<class T>
  Dglap<T>::Dglap(std::function<Set<Operator>(int const&, double const&)> const& SplittingFunctions,
                  std::function<Set<Operator>(bool const&, int const&)>   const& MatchingConditions,
                  Set<T>                                                  const& ObjRef,
                  double                                                  const& MuDistRef,
                  std::vector<double>                                     const& Thresholds,
                  int                                                     const& nsteps):
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
    Set<Operator> MC = _MatchingConditions(Up, (Up ? nf : nf - 1));

    // Convolute matching conditions with a set of objects having the
    // same convolution map of the matching conditions but containing
    // the same objects of the input set of functions 'f'.
    Set<T> MO = MC * Set<T> {MC.GetMap(), f.GetObjects()};

    // Set for 'MO' the convolution map of the next evolution step
    MO.SetMap(_SplittingFunctions((Up ? nf + 1 : nf - 1), 1).GetMap());

    return MO;
  }

  //_________________________________________________________________________________
  template<class T>
  Set<T> Dglap<T>::Derivative(int const& nf, double const& t, Set<T> const& f) const
  {
    return _SplittingFunctions(nf, exp(t / 2)) * f;
  }

  // Fixed template types.
  template class Dglap<Distribution>;
  template class Dglap<Operator>;

  //_________________________________________________________________________________
  template<>
  void Dglap<Distribution>::SetInitialDistributions(std::function<double(int const&, double const&)> const& InDistFunc)
  {
    // Allocate initial scale distributions.
    std::map<int,Distribution> DistMap;
    for (int i = 0; i <= 12; i++)
      DistMap.insert({i, Distribution{_ObjRef.at(0).GetGrid(), InDistFunc, i}});

    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution> {_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistMap});
  }

  //_________________________________________________________________________________
  template<>
  void Dglap<Distribution>::SetInitialDistributions(std::function<std::map<int, double>(double const&)> const& InDistFunc)
  {
    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution> {_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistributionMap(_ObjRef.at(0).GetGrid(), InDistFunc)});
  }

  //_________________________________________________________________________________
  template<>
  void Dglap<Distribution>::SetInitialDistributions(std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc, double const& mu)
  {
    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution> {_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistributionMap(_ObjRef.at(0).GetGrid(), InDistFunc, mu)});
  }
}

