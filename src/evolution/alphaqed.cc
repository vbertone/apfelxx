//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqed.h"
#include "apfel/constants.h"
#include "apfel/betaqed.h"
#include "apfel/messages.h"
#include "apfel/tools.h"

#include <stdexcept>

namespace apfel {
  //_________________________________________________________________________________
  AlphaQED::AlphaQED(double              const& AlphaRef,
		     double              const& MuRef,
		     std::vector<double> const& QuarkThresholds,
		     std::vector<double> const& LeptThresholds,
		     int                 const& pt,
		     int                 const& nstep):
    MatchedEvolution(AlphaRef, MuRef, ConcatenateAndSortVectors(QuarkThresholds, LeptThresholds), nstep),
    _pt(pt)
  {
    // Beta function lambda function.
    _BetaFunction = [&] (int const& nfl, double const& a)-> double
      {
	const double Qr = ConcatenateAndSortVectors(QuarkThresholds, LeptThresholds)[nfl-1] + eps8;
	const int nf = NF(Qr, QuarkThresholds);
	const int nl = NF(Qr, LeptThresholds);
	double bt = 0, powa = a * a;
	for (int i = 0; i <= _pt; i++)
	  {
	    bt -= powa * betaQED(i, nf, nl);
	    powa *= a;
	  }
	return bt;
      };
  }

  //_________________________________________________________________________________
  double AlphaQED::MatchObject(bool const&, int const&, double const& Coup) const
  {
    return Coup;
  }

  //_________________________________________________________________________________
  double AlphaQED::Derivative(int const& nfl, double const&, double const& a) const
  {
    return _BetaFunction(nfl, a);
  }

  //_________________________________________________________________________________
  double AlphaQED::betaQED(int const& pt, int const& nf, int const& nl) const
  {
    double res;
    if (pt == 0)
      res = beta0qed(nf, nl);
    else if (pt == 1)
      res = beta1qed(nf, nl);
    else
      throw std::runtime_error(error("AlphaQED::betaQED","perturbive range out-of-range."));

    return res / pow(FourPi,pt+1);
  }
}
