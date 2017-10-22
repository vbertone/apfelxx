//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/alphaqcd.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/messages.h"
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
    MatchedEvolution(AlphaRef, MuRef, Thresholds, nstep),
    _pt(pt)
  {
    // Initialize all coefficients of the QCD beta function for all
    // numbers of flavours.
    _bQCD.resize(4,3);
    for (auto ipt = 0; ipt <= 2; ipt++)
	for (auto nf = 3; nf <= 6; nf++)
	_bQCD(nf-3, ipt) = betaQCD(ipt, nf);

    // Compute logs of muth2 / m2 needed by the matching conditions.
    vector<double> LogKth;
    for (int im = 0; im < (int) Thresholds.size(); im++)
      if (Thresholds[im] < eps12 || Masses[im] < eps12)
	LogKth.push_back(0);
      else
	LogKth.push_back(2 * log( Thresholds[im] / Masses[im] ));

    // Beta function lambda function.
    _BetaFunction = [&] (int const& nf, double const& as)-> double
      {
	double bt = 0, powas = as * as;
	for (auto i = 0; i <= _pt; i++)
	  {
	    bt -= powas * _bQCD(nf-3, i);
	    powas *= as;
	  }
	return bt;
      };

    // Matching condition lambda function.
    _MatchingConditions = [&,LogKth] (bool const& Up, int const& nf, double const& Coup)-> double
      {
	const auto sgn = ( Up ? 1 : -1);
	const auto ep = Coup / FourPi;
	const double c[3] = { 1, sgn * 2. / 3. * LogKth[nf-1], 4. / 9. * pow(LogKth[nf-1],2) + sgn *  38. / 3. * LogKth[nf-1] + sgn * 14. / 3. };
	double match = 0, powep = 1;
	for (auto i = 0; i <= _pt; i++)
	  {
	    match += c[i] * powep;
	    powep *= ep;
	  }
	return Coup * match;
      };
  }

  //_________________________________________________________________________________
  AlphaQCD::AlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, int const& pt, int const& nstep):
    AlphaQCD(AlphaRef, MuRef, Masses, Masses, pt, nstep)
  {
  }

  //_________________________________________________________________________________
  double AlphaQCD::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    return _MatchingConditions(Up, nf, Coup);
  }

  //_________________________________________________________________________________
  double AlphaQCD::Derivative(int const& nf, double const&, double const& as) const
  {
    return _BetaFunction(nf, as);
  }

  //_________________________________________________________________________________
  double AlphaQCD::betaQCD(int const& pt, int const& nf) const
  {
    double res;
    if ( pt == 0 )
      res = beta0(nf);
    else if ( pt == 1 )
      res = beta1(nf);
    else if ( pt == 2 )
      res = beta2(nf);
    else
      throw runtime_error(error("AlphaQCD::betaQCD","perturbive range out-of-range."));
    return res / pow(FourPi,pt+1);
  }
}
