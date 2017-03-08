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
    MatchedEvolution(ObjRef, MuDistRef, Masses, Thresholds),
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
  Set<Distribution> Dglap::EvolveObject(int const& nf, double const& mu02, double const& mu2, Set<Distribution> const& sd0) const
  {
    // Return immediately "sd0" if "mu02" and "mu2" are equal
    if (mu02 == mu2)
       return sd0;

    // Numerical solution of the evolution equation with fourth-order Runge-Kutta.
    const auto df = rk4<Set<Distribution>>([&](double const& t, Set<Distribution> const& f)->Set<Distribution>{ return Derivative(nf, t, f); });

    // Use "_nstep" steps for the evolution.
    auto f = sd0;
    const auto dt = log( mu2 / mu02 ) / _nstep;
    auto t = log(mu02);
    for (auto k = 0; k < _nstep; k++)
      {
	f += df(t, f, dt);
	t += dt;
      }

    return f;
  }

  //_________________________________________________________________________________
  Set<Distribution> Dglap::MatchObject(bool const& Up, int const& nf, Set<Distribution> const& f) const
  {
    auto MC  = _MatchingConditions(Up, nf, _LogTh2M2[nf]);
    auto MC1 = _MatchingConditions(Up, nf+1, _LogTh2M2[nf]);
    auto MO  = MC * f;
    return Set<Distribution>{MC1.GetMap(), MO.GetObjects()};
  }

  //_________________________________________________________________________________
  Set<Distribution> Dglap::Derivative(int const& nf, double const& t, Set<Distribution> const& f) const
  {
    return _SplittingFunctions(nf, exp(t/2)) * f;
  }

}

