//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <cmath>

// Test expression with singularity at y = 1 to be integrated between
// 0 and infinity and to be treated as a principal value
class PrincipalValueERBL: public apfel::Expression
{
public:
  PrincipalValueERBL(): Expression() {}
  double Singular(double const& y) const
  {
    return 1 / ( 1 - y );
  }
  double LocalPP(double const& y)    const
  {
    return log( 1 - y );
  }
};

int main()
{
  // Grid
  const apfel::Grid g{{{200, 1e-9, 3}, {300, 1e-1, 3}, {1000, 9e-1, 3}}};

  // Distribution
  const apfel::Distribution d = apfel::Operator{g, PrincipalValueERBL{}, apfel::eps5, true} * apfel::Distribution{g, [&] (double const& y) -> double{ return y; }};

  // Tabulation parameters
  const int nx = 100;
  const double xmin = 1e-5;
  const double xmax = 0.999;
  const double xstp = exp( log(xmax / xmin) / ( nx - 1 ) );

  // Check the numerical accuracy of the ERBL principal value
  std::cout << "\nChecking ERBL-like principal value ... " << std::endl;
  std::cout << "    x             numerical       analytic         ratio" << std::endl;
  std::cout << std::scientific;
  for (double x = xmin; x < xmax * 1.000001; x *= xstp)
    {
      const double num = d.Evaluate(x) / x;
      const double ana = log( ( 1 - x ) / x );
      std::cout << x << "\t" << num << "\t" << ana << "\t" << num / ana << std::endl;
    }

  return 0;
}
