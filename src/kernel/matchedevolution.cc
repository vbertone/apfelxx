//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchedevolution.h"
#include "apfel/constants.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/ode.h"
#include "apfel/doubleobject.h"
#include "apfel/doubledistribution.h"
#include "apfel/doubleoperator.h"
#include "apfel/distributionoperator.h"
#include "apfel/operatordistribution.h"

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
    // Compute squared reference scale
    _MuRef2 = pow(MuRef, 2);

    // Compute log of the squared final scale
    _LogMuRef2 = log(_MuRef2);

    // Compute squared thresholds
    for (auto const& th : Thresholds)
      {
        const double th2 = pow(th, 2);
        _Thresholds2.push_back(th2);
        _LogThresholds2.push_back(( th2 > 0 ? log(th2) : -100));
      }

    // Sort the quark thresholds and thair logarithm
    if (_Thresholds2.size() > 1)
      std::sort(_Thresholds2.begin(), _Thresholds2.end());
  }

  //_________________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::EvolveObject(int const& nf, double const& t0, double const& t1, T const& Obj0) const
  {
    // Return immediately "Obj0" if "t0" and "t1" are equal
    if (t0 == t1)
      return Obj0;

    // Numerical solution of the evolution equation with fourth-order
    // Runge-Kutta.
    const auto dObj = rk4<T>([&] (double const& t, T const& Obj) -> T{ return Derivative(nf, t, Obj); });

    // Use "_nsteps" steps for the evolution
    double t   = t0;
    T      Obj = Obj0;
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
    const double mu2  = pow(mu, 2);
    const double lmu2 = log(mu2);

    // Find initial and final number of flavours
    const int nfi = NF(_MuRef2, _Thresholds2);
    const int nff = NF(mu2, _Thresholds2);

    // Do not do the matching if initial and final numbers of flavours
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
    for (int inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
        // Final scale
        const double tf = _LogThresholds2[(sgn ? inf : inf - 1)];

        // Do the matching
        vobj = MatchObject(sgn, inf, EvolveObject(inf, ti, tf, vobj));

        // Update initial scale and displace it by "eps8" to make sure
        // to be above (below) the threshold.
        ti = tf * ( 1 + (sgn ? 1 : -1) * eps8 );
      }
    return EvolveObject(nff, ti, lmu2, vobj);
  }

  // Specialisations
  //_________________________________________________________________________________
  template class MatchedEvolution<double>;
  template class MatchedEvolution<matrix<double>>;
  template class MatchedEvolution<Distribution>;
  template class MatchedEvolution<Set<Distribution>>;
  template class MatchedEvolution<DoubleObject<Distribution>>;
  template class MatchedEvolution<Operator>;
  template class MatchedEvolution<Set<Operator>>;
  template class MatchedEvolution<DoubleObject<Operator>>;
  template class MatchedEvolution<DoubleObject<Distribution, Operator>>;
  template class MatchedEvolution<DoubleObject<Operator, Distribution>>;
  template class MatchedEvolution<Set<DoubleObject<Distribution, Operator>> >;
  template class MatchedEvolution<Set<DoubleObject<Operator, Distribution>> >;
  template class MatchedEvolution<DoubleDistribution>;
  template class MatchedEvolution<DoubleOperator>;
  template class MatchedEvolution<DistributionOperator>;
  template class MatchedEvolution<OperatorDistribution>;
}
