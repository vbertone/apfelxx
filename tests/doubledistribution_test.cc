//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  // Grids
  const apfel::Grid g1{{{100, 5e-6, 3}, {100, 1e-1, 3}, {40, 8e-1, 3}}};
  const apfel::Grid g2{{{100, 5e-6, 3}, {100, 1e-1, 3}, {40, 8e-1, 3}}};

  // Function to interpolate
  const std::function<double(double const&)> fx1 = [] (double const& x1) -> double { return pow(1 - x1, 3); };
  const std::function<double(double const&)> fx2 = [] (double const& x2) -> double { return pow(1 - x2, 2); };
  const std::function<double(double const&, double const&)> fx1x2 = [=] (double const& x1, double const& x2) -> double { return fx1(x1) * fx2(x2); };

  // Derivative
  const std::function<double(double const&, double const&)> dfx1x2 = [=] (double const& x1, double const& x2) -> double { return 6 * pow(1 - x1, 2) * ( 1 - x2 ); };

  // Indefinite integral
  const double xup = 0.01;
  const std::function<double(double const&, double const&)> Sfx1x2 = [=] (double const& x1, double const& x2) -> double
  {
    return -0.08333333333333333*((-pow(-1 + x1,4) + pow(-1 + xup,4))*(-(x2*(3 + (-3 + x2)*x2)) + xup*(3 + (-3 + xup)*xup)));
  };

  // Double distribution
  const apfel::DoubleDistribution dd{g1, g2, fx1x2};

  // Print distribution
  //std::cout << dd;

  // Interpolation test
  const int nx = 10;
  const double xmin = 1e-5;
  const double xmax = 0.95;
  const double xstp = exp(log(xmax / xmin) / ( nx - 1 ));

  std::cout << std::scientific;
  std::cout << "\nChecking the numerical accuracy of the DoubleDistribution object ... " << std::endl;
  std::cout << "        x1      "
            << "        x2      "
            << "  an. func.     "
            << "inter. func.    "
            << "    ratio       "
            << "  an. deriv.    "
            << "inter. deriv.   "
            << "    ratio       "
            << "  an. integ.    "
            << "inter. integ.   "
            << "    ratio       "
            << std::endl;
  for (double x1 = xmin; x1 <= 1.00001 * xmax; x1 *= xstp)
    for (double x2 = xmin; x2 <= 1.00001 * xmax; x2 *= xstp)
      std::cout << x1 << "\t" << x2 << "\t"
                << fx1x2(x1, x2) << "\t" << dd.Evaluate(x1, x2) << "\t" << dd.Evaluate(x1, x2) / fx1x2(x1, x2) << "\t"
                << dfx1x2(x1, x2) << "\t" << dd.Derive(x1, x2) << "\t" << dd.Derive(x1, x2) / dfx1x2(x1, x2) << "\t"
                << Sfx1x2(x1, x2) << "\t" << dd.Integrate(x1, xup, x2, xup) << "\t" << dd.Integrate(x1, xup, x2, xup) / Sfx1x2(x1, x2) << "\t"
                << std::endl;

  // Single distribution
  const apfel::Distribution d1{g1, fx1};
  const apfel::Distribution d2{g2, fx2};

  // Double distribution
  const apfel::DoubleDistribution dds{d1, d2};

  std::cout << "\nChecking the numerical accuracy of the DoubleDistribution object constructed using two single distributions ... " << std::endl;
  std::cout << "        x1      "
            << "        x2      "
            << "  an. func.     "
            << "inter. func.    "
            << "    ratio       "
            << std::endl;
  for (double x1 = xmin; x1 <= 1.000001 * xmax; x1 *= xstp)
    for (double x2 = xmin; x2 <= 1.000001 * xmax; x2 *= xstp)
      std::cout << x1 << "\t" << x2 << "\t" << fx1(x1) * fx2(x2) << "\t" << dds.Evaluate(x1, x2) << "\t" << dds.Evaluate(x1, x2) / ( fx1(x1) * fx2(x2) ) << std::endl;
  std::cout << "\n";

  apfel::Timer t;
  const int nint = 1000000;
  std::cout << "Performance test ("<< nint << " interpolations)... ";
  for (int i = 0; i < nint; i++)
    dd.Evaluate(0.1111, 0.2222);
  t.stop();
  t.start();
  std::cout << "Performance test ("<< nint << " derivations)... ";
  for (int i = 0; i < nint; i++)
    dd.Evaluate(0.1111, 0.2222);
  t.stop();
  t.start();
  std::cout << "Performance test ("<< nint / 1000 << " integrations)... ";
  for (int i = 0; i < nint / 1000; i++)
    dd.Integrate(0.1111, 0.55555, 0.1111, 0.55555);
  t.stop();
  std::cout << "\n";

  return 0;
}
