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
    MatchedEvolution(AlphaRef, MuRef, Masses, Thresholds, nstep),
    _pt(pt)
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
  double AlphaQCD::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    const auto sgn = ( Up ? 1 : -1);
    const auto ep = Coup / FourPi;
    const auto LogKth = _LogTh2M2[nf];
    const double c[3] = { 1, sgn * 2. / 3. * LogKth, 4. / 9. * pow(LogKth,2) + sgn *  38. / 3. * LogKth + sgn * 14. / 3. };
    double match = 0, powep = 1;
    for (auto i = 0; i <= _pt; i++)
      {
        match += c[i] * powep;
        powep *= ep;
      }
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
    return res / pow(FourPi,pt+1);
  }

  //_________________________________________________________________________________
  double AlphaQCD::Derivative(int const& nf, double const&, double const& as) const
  {
    double bt = 0, powas = as * as;
    for (auto i = 0; i <= _pt; i++)
      {
        bt -= powas * _bQCD[nf-3][i];
        powas *= as;
      }
    return bt;
  }

}
