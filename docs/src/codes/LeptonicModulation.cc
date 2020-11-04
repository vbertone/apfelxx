/*
 * Author: Valerio Bertone
 */

// APFEL++ libs
#include "apfel/apfelxx.h"

int main() {
  // Define vector-boson kinematics
  const double Q  = 91;
  const double y  = 0;
  const double Q2 = Q * Q;
  const std::vector<double> qTv = {1, 2, 3, 4, 5};

  for (auto const& qT : qTv)
    {
      const double qT2 = qT * qT;
      const double M2  = Q2 + qT2;
      const double M   = sqrt(M2);
      const double b   = qT / M;
      const double b2  = b * b;

      const int npT = 10000;
      const double pTmin = 41.5;
      const double pTmax = 45.5;
      const double step  = ( pTmax - pTmin ) / ( npT - 1 );

      double pT = pTmin;
      std::cout << std::scientific;
      for (int ipT = 0; ipT < npT; ipT++)
	{
	  // Define auxiliary variables
	  const double a  = Q / 2 / pT;
	  const double c  = Q2 / 2 / pT / M;
	  const double a2 = a * a;

	  // Define integrand
	  const auto Integrand = [=] (double const& y) -> double
	    {
	      if ((y < 1 + c && y > - 1 + c) || (y < 1 - c && y > - 1 - c))
		return 0;

	      const double fp = ( ( a2 - 1 ) + pow(y + c, 2) ) / sqrt( pow(y + c, 2) - 1 );
	      const double fm = ( ( a2 - 1 ) + pow(y - c, 2) ) / sqrt( pow(y - c, 2) - 1 );
	      return ( fp + fm ) / sqrt(b2 - pow(y, 2));
	    };

	  // Define integral
	  const apfel::Integrator Integral{Integrand};

	  // Integrate and print result
	  std::cout << qT << "\t" << pT << "\t" << 3 * pT * Integral.integrate(-b, b, 5e-5) / 2 / M_PI / Q2 << std::endl;

	  // Increment step
	  pT += step;
	}
      std::cout << "\n\n";
    }
  return 0;
}
