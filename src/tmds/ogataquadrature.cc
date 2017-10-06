//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/ogataquadrature.h"

using namespace std;

namespace apfel {
  //_____________________________________________________________________________
  double OgataQuadrature(function<double(double const&)> const& func, double const& qT, double const& CutOff, double const& h)
  {
    // Define helper functions.
    const auto psi  = [] (double const& t)->double{ return t * tanh( M_PI * sinh(t) / 2 ); };
    const auto psip = [] (double const& t)->double{ return ( M_PI * t * cosh(t) + sinh( M_PI * sinh(t) ) ) / ( 1 + cosh( M_PI * sinh(t) ) ); };

    // Number of terms counter.
    int i = 0;

    // Run over the zero's of J0. There are 1000 zeros
    // percalculated. This should be enough for all practical
    // applications.
    double integral = 0;
    for (auto const& jz : j0Zeros)
      {
	i++;
	const double w = y0(jz) / j1(jz);
	const double z = jz / M_PI;
	const double x = M_PI * psi( h * z ) / h;
	const double f = x * func( x / qT ) / 2;
	const double term = f * j0(x) * psip( h * z );

	// Break when the absolute value of the last term is less
	// than "CutOff". This assumes that terms are increasingly
	// small. Provided that the integrand is not badly behaved,
	// this should be guaranteed by the Bessel function j0.
	if (abs(term) < CutOff)
	  break;

	integral += w * term;
      }
    integral *= M_PI / qT / qT;

    // If the number of terms is equal to the size of the "j0Zeros"
    // vector give a worning.
    if (i == (int) j0Zeros.size())
      warning("OgataQuadrature", "Number of j0 zero's available exceeded: the integration might not have converged.");

    return integral;
  }
}
