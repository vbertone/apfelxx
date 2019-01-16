/*
  Author: Valerio Bertone
 */

// Standard libs

// APFEL++ libs
#include "apfel/apfelxx.h"

// b* prescription
double bstar(double const& b)
{
  const double bmax = 2 * exp( - apfel::emc);
  return b / sqrt( 1 + pow(b / bmax, 2) );
}

// Non-perturnative function
double fNP(double const& b, double const& zetaf)
{
  const double g1 = 0.02;
  const double g2 = 0.5;
  const double Q0 = 1.6;
  return exp( ( - g1 - g2 * log( sqrt(zetaf) / Q0 / 2 ) ) * b * b / 2 );
}

int main() {
  // Ogata-quadrature object of degree zero.
  apfel::OgataQuadrature OgataObj{0, 1e-11};

  const double Q = 10;

  for (int ord = 0; ord < 4; ord++)
    {
      const auto integrand = [=] (double const& b) -> double
	{
	  const double f1 = fNP(b, Q * Q);
	  const double f2 = f1;
	  const double lg = 2 * log(bstar(b) * Q / 2 / exp( - apfel::emc));
	  return b * f1 * f2 * pow(lg, ord) / 2;
	};

      const int nqT = 1000;
      const double qTmin = 1e-2;
      const double qTmax = 5 * Q;
      const double qTstep = exp( log( qTmax / qTmin ) / ( nqT - 1 ) );

      double qT = qTmin;
      for (int iqT = 0; iqT < nqT; iqT++)
	{
	  // Perform Fourier transform and obtain cross section
	  const double trans = OgataObj.transform(integrand, qT);
	  double exact = 0;
	  if (ord == 1)
	    exact = - 1 / qT / qT;
	  else if (ord == 2)
	    exact = - 4 * log(Q / qT) / qT / qT;
	  else if (ord == 3)
	    exact = - 12 * pow(log(Q / qT), 2) / qT / qT;

	  const double powsup  = pow(qT/Q, 0.2);
	  const double matched = ( 1 - powsup ) * exact + powsup * trans;

	  // Print results
	  std::cout << std::scientific << ord << "\t" << qT << "\t" << trans << "\t" << exact << "\t" << matched << std::endl;

	  qT *= qTstep;
	}
      std::cout << std::endl;
    }

  return 0;
}
