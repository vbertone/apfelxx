//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <algorithm>
#include <cmath>

#include "apfel/coupling.h"
#include "apfel/tools.h"

using namespace std;

namespace apfel
{
  //_________________________________________________________________________
  Coupling::Coupling(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, vector<double> const& Thresholds):
    _AlphaRef(AlphaRef)
  {
    // Compute squared final scale
    _MuRef2 = pow(MuRef,2);

    // Compute squared thresholds
    for (auto &th : Thresholds)
      _Thresholds2.push_back(pow(th,2));

    // Check that "Masses" and "Thresholds" have the same size
    if (Masses.size() != Thresholds.size())
      throw logic_exception("Coupling::Coupling", "Masses and Thresholds vectors have diffrent sizes.");

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
  Coupling::Coupling(double const& AlphaRef, double const& MuRef, vector<double> const& Masses):
    Coupling(AlphaRef, MuRef, Masses, Masses)
  {
  }

  //_________________________________________________________________________
  double Coupling::GetCoupling(double const& mu) const
  {
    auto const mu2 = pow(mu,2);

    // Find initial and final number of flavours
    const auto nfi = lower_bound(_Thresholds2.begin()+1, _Thresholds2.end(), _MuRef2) - _Thresholds2.begin();
    const auto nff = lower_bound(_Thresholds2.begin()+1, _Thresholds2.end(), mu2)     - _Thresholds2.begin();
    const auto sgn = signbit(nfi - nff);

    if ( nfi == nff )
      return Coup(nfi, _AlphaRef, _MuRef2, mu2);

    auto asi  = _AlphaRef;
    auto asf  = _AlphaRef;
    auto mui2 = _MuRef2;
    auto muf2 = _Thresholds2[nfi];
    for(auto inf = nfi; (sgn ? inf < nff : inf > nff); inf += (sgn ? 1 : -1))
      {
	asf  = Coup(inf, asi, mui2, muf2);
	asf  = MatchCoupling(sgn, asf, _LogTh2M2[inf]);
	asi  = asf;
	mui2 = muf2;
	muf2 = _Thresholds2[fmin(inf+1,nff)];
      }
    return Coup(nff, asf, mui2, mu2);
  }

}
