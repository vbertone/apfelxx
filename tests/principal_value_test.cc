//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

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
  double LocalPP(double const& y) const
  {
    return log( 1 - y );
  }
};

// Test expression with singularity at y = 1/kappa to be integrated
// between 0 and 1 and to be treated as a principal value
class PrincipalValueDGLAP: public apfel::Expression
{
public:
  PrincipalValueDGLAP(double const& xi): Expression(1/xi) {}
  double SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 1 / ( 1 - kappa * y );
    else
      return 0;
  }
  double LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return log( 1 - kappa * y ) / kappa;
    else
      return 0;
  }
};

int main()
{
  // Grid
  const apfel::Grid g{{{1000, 1e-7, 3}, {1000, 1e-1, 3}, {1000, 9e-1, 3}}};

  // Define skewness
  const double xi = 0.09;

  // Distribution for principa-value at y = 1 integrated between x and
  // infinity
  const apfel::Distribution d1 = apfel::Operator{g, PrincipalValueERBL{}, apfel::eps5, true} * apfel::Distribution{g, [&] (double const& y) -> double{ return y; }};

  // Distribution for principa-value at y = 1/kappa integrated between
  // x and 1
  const apfel::Distribution dk = apfel::Operator{g, PrincipalValueDGLAP{xi}, apfel::eps5, true} * apfel::Distribution{g, [&] (double const&) -> double{ return 1; }};

  // Tabulation parameters
  const int nx = 100;
  const double xmin = 1e-5;
  const double xmax = 0.99;
  const double xstp = exp( log(xmax / xmin) / ( nx - 1 ) );

  // Check the numerical accuracy of the ERBL principal value
  std::cout << "\nChecking ERBL-like principal value ... " << std::endl;
  std::cout << "    x             numerical       analytic         ratio" << std::endl;
  std::cout << std::scientific;
  for (double x = xmin; x < xmax * 1.000001; x *= xstp)
    {
      const double num = d1.Evaluate(x) / x;
      const double ana = log( ( 1 - x ) / x );
      std::cout << x << "\t" << num << "\t" << ana << "\t" << num / ana << std::endl;
    }

  // Check the numerical accuracy of the DGLAP principal value
  std::cout << "\nChecking DGLAP-like principal value ... " << std::endl;
  std::cout << "    x             numerical       analytic         ratio" << std::endl;
  std::cout << std::scientific;
  for (double x = xmin; x < xmax * 1.000001; x *= xstp)
    {
      const double kappa = xi / x;
      const double num = dk.Evaluate(x);
      const double ana = (xi > x ? x / xi * log( x * ( 1 - xi ) / ( xi - x ) ): 0);
      const apfel::Integrator I{[&] (double const& y) -> double { return ( kappa * y > 1 ? - 1 / kappa / y : 0 ); }};
      std::cout << x << "\t"
                << num << "\t" << ana << "\t" << num / ana << "\t"
                << ( I.integrate(x, 1, apfel::eps9) + log( kappa * ( 1 - kappa * x ) / ( kappa - 1 ) ) / kappa ) / ana
                << std::endl;
    }

  return 0;
}
