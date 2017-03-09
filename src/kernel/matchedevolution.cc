//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/matchedevolution.h"
#include "apfel/tools.h"
#include "apfel/distribution.h"
#include "apfel/set.h"
#include "apfel/ode.h"

#include <algorithm>

using std::min;

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  MatchedEvolution<T>::MatchedEvolution(T              const& ObjRef,
					double         const& MuRef,
					vector<double> const& Masses,
					vector<double> const& Thresholds,
					int            const& nsteps):
    _ObjRef(ObjRef),
    _Masses(Masses),
    _Thresholds(Thresholds),
    _nsteps(nsteps)
  {
    // Check that "Masses" and "Thresholds" have the same size
    if (Masses.size() != Thresholds.size())
      throw logic_exception("MatchedEvolution::MatchedEvolution", "Masses and Thresholds vectors have diffrent sizes.");

    // Compute squared reference scale
    _MuRef2 = pow(MuRef,2);

    // Compute log of the squared final scale
    _LogMuRef2 = log(_MuRef2);

    // Compute squared thresholds
    for (auto &th : Thresholds)
      {
	const double th2 = pow(th,2);
	_Thresholds2.push_back(th2);
	_LogThresholds2.push_back(( th2 > 0 ? log(th2) : -100));
      }

    // Compute logs of muth2 / m2
    for (auto im = 0; im < (int) Thresholds.size(); im++)
      if (_Thresholds2[im] == 0 || Masses[im] == 0)
	_LogTh2M2.push_back(-100);
      else
	_LogTh2M2.push_back(log(_Thresholds2[im] / pow(Masses[im],2)));

    // Sort the quark thresholds and logs
    if (_Thresholds2.size() > 1)
      sort(_Thresholds2.begin(), _Thresholds2.end());
  }

  //_________________________________________________________________________
  template<class T>
  MatchedEvolution<T>::MatchedEvolution(T              const& ObjRef,
					double         const& MuRef,
					vector<double> const& Masses,
					int            const& nsteps):
    MatchedEvolution(ObjRef, MuRef, Masses, Masses, nsteps)
  {
  }


  //_________________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::EvolveObject(int const& nf, double const& t0, double const& t1, T const& Obj0) const
  {
    // Return immediately "Obj0" if "t0" and "t1" are equal
    if (t0 == t1)
      return Obj0;

    // Numerical solution of the evolution equation with fourth-order Runge-Kutta.
    const auto dObj = rk4<T>([&](double const& t, T const& Obj)->T{ return Derivative(nf, t, Obj); });

    // Use "_nsteps" steps for the evolution.
    auto t        = t0;
    auto Obj      = Obj0;
    const auto dt = ( t1 - t0 ) / _nsteps;
    for (auto k = 0; k < _nsteps; k++)
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
    auto const mu2  = pow(mu,2);
    auto const lmu2 = log(mu2);

    // Find initial and final number of flavours
    const auto nfi = lower_bound(_Thresholds2.begin()+1, _Thresholds2.end(), _MuRef2) - _Thresholds2.begin();
    const auto nff = lower_bound(_Thresholds2.begin()+1, _Thresholds2.end(),     mu2) - _Thresholds2.begin();

    // Don't do the matching is initial and final number of flavours are equal
    if ( nfi == nff )
      return EvolveObject(nfi, _LogMuRef2, lmu2, _ObjRef);

    // Direction of the evolution
    const auto sgn = std::signbit(nfi - nff);

    // Create a vector of objects containing the object right above each threshold
    // to make sure that every time a threshold is crossed a new object with a
    // different convolution map is created (effective only when a "Set" object
    // is evolved).
    vector<T> vobj = { _ObjRef };
    auto ti        = _LogMuRef2;
    auto tf        = _LogThresholds2[nfi];
    for(auto inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
	vobj.push_back(MatchObject(sgn, inf, EvolveObject(inf, ti, tf, vobj.back())));
	ti = tf + eps8;                    // Add "eps8" to make sure to be above the threshold
	tf = _LogThresholds2[min(inf+1,nff-1)];
      }
    return EvolveObject(nff, ti, lmu2, vobj.back());
  }

  // template fixed types
  template class MatchedEvolution<double>;               //<! Single coupling
  template class MatchedEvolution<Set<Distribution>>;    //<! Set of distributions

}
