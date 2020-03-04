//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchedevolution.h"
#include "apfel/constants.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/ode.h"
#include "apfel/doubleobject.h"

#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  MatchedEvolution<T>::MatchedEvolution(T                   const& ObjRef,
                                        double              const& MuRef,
                                        std::vector<double> const& Thresholds,
                                        int                 const& nsteps):
    _ObjRef(ObjRef),
    _MuRef(MuRef),
    _Thresholds(Thresholds),
    _nsteps(nsteps)
  {
    // Compute squared reference scale.
    _MuRef2 = pow(MuRef,2);

    // Compute log of the squared final scale.
    _LogMuRef2 = log(_MuRef2);

    // Compute squared thresholds.
    for (auto const& th : Thresholds)
      {
        const double th2 = pow(th, 2);
        _Thresholds2.push_back(th2);
        _LogThresholds2.push_back(( th2 > 0 ? log(th2) : -100));
      }

    // Sort the quark thresholds and logs.
    if (_Thresholds2.size() > 1)
      std::sort(_Thresholds2.begin(), _Thresholds2.end());
  }

  //_________________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::EvolveObject(int const& nf, double const& t0, double const& t1, T const& Obj0) const
  {
    // Return immediately "Obj0" if "t0" and "t1" are equal.
    if (t0 == t1)
      return Obj0;

    // Numerical solution of the evolution equation with fourth-order
    // Runge-Kutta.
    const auto dObj = rk4<T>([&](double const& t, T const& Obj)->T{ return Derivative(nf, t, Obj); });

    // Use "_nsteps" steps for the evolution.
    double t        = t0;
    T      Obj      = Obj0;
    const double dt = ( t1 - t0 ) / _nsteps;
    for (int k = 0; k < _nsteps; k++)
      {
        Obj += dObj(t, Obj, dt);
        t   += dt;
      }
    return Obj;
  }

  //_________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::Evaluate(double const& mu) const
  {
    const double mu2  = pow(mu,2);
    const double lmu2 = log(mu2);

    // Find initial and final number of flavours.
    const int nfi = NF(_MuRef2, _Thresholds2);
    const int nff = NF(mu2, _Thresholds2);

    // Don't do the matching is initial and final number of flavours
    // are equal.
    if (nfi == nff)
      return EvolveObject(nfi, _LogMuRef2, lmu2, _ObjRef);

    // Direction of the evolution
    const bool sgn = std::signbit(nfi - nff);

    // Create a vector of objects containing the object right above
    // each threshold to make sure that every time a threshold is
    // crossed a new object with a different convolution map is
    // created (effective only when a "Set" object is evolved).
    T      vobj = _ObjRef;
    double ti   = _LogMuRef2;
    double tf   = _LogThresholds2[(sgn ? nfi : nfi-1)];
    for (int inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
        vobj = MatchObject(sgn, inf, EvolveObject(inf, ti, tf, vobj));
        ti   = tf + (sgn ? 1 : -1) * eps8;  // Add "eps8" to make sure to be above the threshold
        tf   = (sgn ? _LogThresholds2[std::min(inf+1,nff-1)] : _LogThresholds2[std::max(inf-2,nff-1)]);
      }
    return EvolveObject(nff, ti, lmu2, vobj);
  }

  // template fixed types
  template class MatchedEvolution<double>;                     //<! Single coupling
  template class MatchedEvolution<Distribution>;               //<! Single distribution
  template class MatchedEvolution<Set<Distribution>>;          //<! Set of distributions
  template class MatchedEvolution<DoubleObject<Distribution>>; //<! Set of distributions
  template class MatchedEvolution<Operator>;                   //<! Single Operator
  template class MatchedEvolution<Set<Operator>>;              //<! Set of Operators
}
