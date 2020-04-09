//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <iomanip>

int main()
{
  std::cout << std::setprecision(12) << std::scientific;

  // Grid
  const apfel::Grid g{{{100, 9.9e-6, 5}, {100,1e-1, 5}, {40,8e-1,5}}};

  // Test distribution
  const auto xg = [&] (double const& x) -> double{ return x * ( 1 - x ); };
  const apfel::Distribution xgluon{g, xg};

  // Derivative of the distribution
  const auto dxg = [&] (double const& x) -> double{ return 1 - 2 * x; };

  // Test values
  std::vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.51, 0.6, 0.7, 0.8, 0.9};

  std::cout << "\n        x               "
            << "analytic function       "
            << " inter. function        "
            << "      ratio             "
            << "analytic derivative     "
            << " inter. derivative      "
            << "      ratio             "
            << std::endl;

  for (auto const& ix: x)
    {
      const double original = xg(ix);
      const double derorig  = dxg(ix);
      const double interpol = xgluon.Evaluate(ix);
      const double derive   = xgluon.Derive(ix);
      std::cout << ix << "\t"
                << original << "\t"
                << interpol << "\t"
                << interpol / original<< "\t"
                << derorig  << "\t"
                << derive   << "\t"
                << derive / derorig << std::endl;
    }
  std::cout << "\n";

  apfel::Timer t;
  const int nint = 1000000;
  std::cout << "Performance test ("<< nint << " interpolations)... ";
  for (int i = 0; i < nint; i++)
    xgluon.Evaluate(0.1111);
  t.stop();
  std::cout << "\n";

  return 0;
}
