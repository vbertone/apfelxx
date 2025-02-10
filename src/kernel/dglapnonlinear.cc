//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/dglapnonlinear.h"
#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________________
  DglapNonLinear::DglapNonLinear(std::function<Set<Operator>(int const&, double const&)>                           const& SplittingFunctions,
                                 std::function<Set<Operator>(bool const&, int const&)>                             const& MatchingConditions,
                                 std::function<std::map<int, std::function<double(double const&)>>(double const&)> const& TranformationFuncs,
                                 Set<Distribution>                                                                 const& ObjRef,
                                 double                                                                            const& MuDistRef,
                                 std::vector<double>                                                               const& Thresholds,
                                 int                                                                               const& nsteps):
    MatchedEvolution<Set<Distribution>>(ObjRef, MuDistRef, Thresholds, nsteps),
    _SplittingFunctions(SplittingFunctions),
    _MatchingConditions(MatchingConditions),
    _TranformationFuncs(TranformationFuncs)
  {
  }

  //_________________________________________________________________________________
  Set<Distribution> DglapNonLinear::MatchObject(bool const& Up, int const& nf, Set<Distribution> const& f) const
  {
    // Get matching conditions
    Set<Operator> MC = _MatchingConditions(Up, (Up ? nf : nf - 1));

    // Convolute matching conditions with a set of objects having the
    // same convolution map of the matching conditions but containing
    // the same objects of the input set of functions 'f'.
    Set<Distribution> MO = MC * Set<Distribution> {MC.GetMap(), f.GetObjects()};

    // Set for 'MO' the convolution map of the next evolution step
    MO.SetMap(_SplittingFunctions((Up ? nf + 1 : nf - 1), 1).GetMap());

    return MO;
  }

  //_________________________________________________________________________________
  Set<Distribution> DglapNonLinear::Derivative(int const& nf, double const& t, Set<Distribution> const& f) const
  {
    return _SplittingFunctions(nf, t) * f.Transform(_TranformationFuncs(t));
  }

  //_________________________________________________________________________________
  void DglapNonLinear::SetInitialDistributions(std::function<double(int const&, double const&)> const& InDistFunc)
  {
    // Allocate initial scale distributions
    std::map<int,Distribution> DistMap;
    for (int i = 0; i <= 12; i++)
      DistMap.insert({i, Distribution{_ObjRef.at(0).GetGrid(), InDistFunc, i}});

    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution> {_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistMap});
  }

  //_________________________________________________________________________________
  void DglapNonLinear::SetInitialDistributions(std::function<std::map<int, double>(double const&)> const& InDistFunc)
  {
    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution> {_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistributionMap(_ObjRef.at(0).GetGrid(), InDistFunc)});
  }

  //_________________________________________________________________________________
  void DglapNonLinear::SetInitialDistributions(std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc, double const& mu)
  {
    // Create set of initial distributions (assumed to be in the QCD
    // evolution basis).
    SetObjectRef(Set<Distribution> {_SplittingFunctions(NF(_MuRef, _Thresholds), 0).GetMap(), DistributionMap(_ObjRef.at(0).GetGrid(), InDistFunc, mu)});
  }
}
