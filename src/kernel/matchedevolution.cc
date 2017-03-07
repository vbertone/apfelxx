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

#include <algorithm>
#include <cmath>

using std::min;

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  MatchedEvolution<T>::MatchedEvolution(T const& ObjRef, double const& MuRef, vector<double> const& Masses, vector<double> const& Thresholds):
    _ObjRef(ObjRef),
    _Masses(Masses),
    _Thresholds(Thresholds)
  {
    // Check that "Masses" and "Thresholds" have the same size
    if (Masses.size() != Thresholds.size())
      throw logic_exception("MatchedEvolution::MatchedEvolution", "Masses and Thresholds vectors have diffrent sizes.");

    // Compute squared final scale
    _MuRef2 = pow(MuRef,2);

    // Compute squared thresholds
    for (auto &th : Thresholds)
      _Thresholds2.push_back(pow(th,2));

    // Compute logs of muth2 / m2
    for (auto im = 0; im < (int) Thresholds.size(); im++)
      if (_Thresholds2[im] == 0 || Masses[im] == 0)
	_LogTh2M2.push_back(0);
      else
	_LogTh2M2.push_back(log(_Thresholds2[im] / pow(Masses[im],2)));

    // Sort the quark thresholds and logs
    if (_Thresholds2.size() > 1)
      sort(_Thresholds2.begin(), _Thresholds2.end());
  }

  //_________________________________________________________________________
  template<class T>
  MatchedEvolution<T>::MatchedEvolution(T const& ObjRef, double const& MuRef, vector<double> const& Masses):
    MatchedEvolution(ObjRef, MuRef, Masses, Masses)
  {
  }

  //_________________________________________________________________________
  template<class T>
  T MatchedEvolution<T>::GetObject(double const& mu) const
  {
    auto const mu2 = pow(mu,2);

    // Find initial and final number of flavours
    const auto nfi = lower_bound(_Thresholds2.begin()+1, _Thresholds2.end(), _MuRef2) - _Thresholds2.begin();
    const auto nff = lower_bound(_Thresholds2.begin()+1, _Thresholds2.end(), mu2)     - _Thresholds2.begin();

    // Don't do the matching is initial and final number of flavours are equal
    if ( nfi == nff )
      return EvolveObject(nfi, _MuRef2, mu2, _ObjRef);

    // Direction of the evolution
    const auto sgn = std::signbit(nfi - nff);

    // Create a vector of objects containing the object right above each threshold
    // to make sure that every time a threshold is crossed a new object with a
    // different convolution map is created (effective only when a "Set" object
    // is evolved).
    vector<T> vobj = { _ObjRef };
    auto mui2      = _MuRef2;
    auto muf2      = _Thresholds2[nfi];
    for(auto inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
	vobj.push_back(MatchObject(sgn, inf, EvolveObject(inf, mui2, muf2, vobj.back())));
	mui2 = muf2 + eps12;
	muf2 = _Thresholds2[min(inf+1,nff-1)];
      }
    return EvolveObject(nff, mui2, mu2, vobj.back());
  }

  // template fixed types
  template class MatchedEvolution<double>;               //<! Single coupling
  template class MatchedEvolution<Distribution>;         //<! Single distribution
  template class MatchedEvolution<Set<Distribution>>;    //<! Set of distributions

}
