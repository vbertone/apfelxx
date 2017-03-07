//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/alphaqcd.h"
#include "apfel/tools.h"
#include "apfel/ode.h"

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  AlphaQCD::AlphaQCD(double         const& AlphaRef,
		     double         const& MuRef,
		     vector<double> const& Masses,
		     vector<double> const& Thresholds,
		     int            const& pt,
		     int            const& nstep):
    MatchedEvolution(AlphaRef, MuRef, Masses, Thresholds),
    _pt(pt),
    _nstep(nstep)
  {
    // Initialize all coefficients of the QCD beta function for all numbers of flavours
    for (auto ipt = 0; ipt <= 2; ipt++)
	for (auto nf = 3; nf <= 6; nf++)
	_bQCD[nf-3][ipt] = betaQCD(ipt, nf);
  }

  //_________________________________________________________________________________
  AlphaQCD::AlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, int const& pt, int const& nstep):
    AlphaQCD(AlphaRef, MuRef, Masses, Masses, pt, nstep)
  {
  }

  //_________________________________________________________________________________
  double AlphaQCD::EvolveObject(int const& nf, double const& mu02, double const& mu2, double const& as0) const
  {
    // Return immediately "as0" if "mu02" and "mu2" are equal
    if (mu02 == mu2)
      return as0;

    // Numerical solution of the evolution equation with fourth-order Runge-Kutta.
    const auto das = rk4([&](double const&, double const& as)->double{ return fbeta(as, _bQCD[nf-3]); });

    // Use "_nstep" steps for the evolution.
    auto as        = as0 / FourPi;
    const auto dlr = log( mu2 / mu02 ) / _nstep;
    for (auto k = 0; k < _nstep; k++)
      as += das(0, as, dlr);

    return FourPi * as;
  }

  //_________________________________________________________________________________
  double AlphaQCD::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    const auto sgn = ( Up ? 1 : -1);
    const auto ep = Coup / FourPi;
    const auto LogKth = _LogTh2M2[nf];
    const double c[] = { 1, sgn * 2. / 3. * LogKth, 4. / 9. * pow(LogKth,2) + sgn *  38. / 3. * LogKth + sgn * 14. / 3. };
    double match = 0;
    for (auto i = 0; i <= _pt; i++) match += c[i] * pow(ep,i);
    return Coup * match;
  }

  //_________________________________________________________________________________
  double AlphaQCD::betaQCD(int const& pt, int const& nf) const
  {
    double res = 0;
    if      ( pt == 0 ) res = ( 33. - 2. * nf ) / 3.;
    else if ( pt == 1 ) res = 102. - 38. / 3. * nf;
    else if ( pt == 2 ) res = 2857. / 2. - 5033. / 18. * nf + 325. / 54. * nf * nf;
    else throw runtime_exception("AlphaQCD::betaQCD","perturbive range out-of-range.");
    return res;
  }

  //_________________________________________________________________________________
  double AlphaQCD::fbeta(double const& as, array<double,3> const& bQCD) const
  {
    double bt = 0, powas = as * as;
    for (auto i = 0; i <= _pt; i++)
      {
        bt -= powas * bQCD[i];
        powas *= as;
      }
    return bt;
  }

}
