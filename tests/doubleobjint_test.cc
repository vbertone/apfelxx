//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  // x-space grid
  //const apfel::Grid g{{{1000, 1e-5, 3}, {200, 1e-1, 3}, {50, 6e-1, 3}, {50, 8e-1, 3}}};
  const apfel::Grid g{{{100, 1e-5, 3}, {60, 1e-1, 3}, {50, 6e-1, 3}, {50, 8e-1, 3}}};

  // This function comes from the Jacobian
  const std::function<double(double const&)> f = [] (double const& x) -> double { return sqrt(x); };

  // Double distribution to be integrated
  apfel::DoubleObject<apfel::Distribution> dd{{{1, apfel::Distribution{g, f}, apfel::Distribution{g, f}}}};

  // Parameters
  const double Vs = 7000;
  const double Q1 = 66;
  const double Q2 = 116;
  const double y1 = 0;
  const double y2 = 2.4;

  // Functions defining the upper and lower integration bound in x2 as
  // functions of x1.
  const double x1min = Q1 / Vs * exp(y1);
  const double x1max = Q2 / Vs * exp(y2);
  const std::function<double(double const&)> x2minf = [=] (double const& x1) -> double { return std::max(pow(Q1 / Vs, 2) / x1, exp(- 2 * y2) * x1); };
  const std::function<double(double const&)> x2maxf = [=] (double const& x1) -> double { return std::min(pow(Q2 / Vs, 2) / x1, exp(- 2 * y1) * x1); };

  // Functions defining the upper and lower integration bound in x1 as
  // functions of x2.
  const double x2min = Q1 / Vs * exp(-y2);
  const double x2max = Q2 / Vs * exp(-y1);
  const std::function<double(double const&)> x1minf = [=] (double const& x2) -> double { return std::max(pow(Q1 / Vs, 2) / x2, exp(2 * y1) * x2); };
  const std::function<double(double const&)> x1maxf = [=] (double const& x2) -> double { return std::min(pow(Q2 / Vs, 2) / x2, exp(2 * y2) * x2); };

  // Print results
  const double an  = ( pow(Q2, 3) - pow(Q1, 3) ) * ( y2 - y1 ) / 3 / pow(Vs, 2);
  const double nm1 = Vs * dd.Integrate(x1min, x1max, x2minf, x2maxf) / 2;
  const double nm2 = Vs * dd.Integrate(x1minf, x1maxf, x2min, x2max) / 2;
  std::cout << "  Analytic       Numerical1      Numerical2        Ratio1          Ratio2" << std::endl;
  std::cout << std::scientific << an << "\t" << nm1 << "\t" << nm2 << "\t" << nm1 / an << "\t" << nm2 / an << std::endl;

  return 0;
}
