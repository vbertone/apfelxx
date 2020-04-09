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
  std::cout << std::setprecision(8) << std::scientific;

  // Test distribution
  const auto xg = [&] (double const& x) -> double{ return x * ( 1 - x ); };

  // Derivative of the distribution
  const auto dxg = [&] (double const& x) -> double{ return 1 - 2 * x; };

  // Integral of the distribution
  const auto ixg = [&] (double const& x) -> double{ return x * x * ( 1. / 2 - x / 3 ); };

  // Grid
  const apfel::Grid g{{{100, 9.9e-6, 5}, {100, 1e-1, 5}, {40, 8e-1, 5}}};

  // Interpolated distribution
  const apfel::Distribution xgluon{g, xg};

  // Test values
  std::vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.51, 0.6, 0.7, 0.8, 0.9};
  //std::vector<double> x = {1e-4};

  std::cout << "\n        x       "
            << "  an. func.     "
            << " inter. func.   "
            << "    ratio       "
            << "  an. deriv.    "
            << " inter. deriv.  "
            << "    ratio       "
            << "  an. integ.    "
            << " inter. integ.  "
            << "    ratio       "
            << std::endl;

  for (auto const& ix: x)
    {
      const double original = xg(ix);
      const double derorig  = dxg(ix);
      const double intorig  = ixg(0.9) - ixg(ix);
      const double interpol = xgluon.Evaluate(ix);
      const double derive   = xgluon.Derive(ix);
      const double integr   = xgluon.Integrate(ix, 0.9);
      std::cout << ix << "\t"
                << original << "\t"
                << interpol << "\t"
                << interpol / original << "\t"
                << derorig  << "\t"
                << derive   << "\t"
                << derive / derorig << "\t"
                << intorig  << "\t"
                << integr   << "\t"
                << integr / intorig << "\t"
                << std::endl;
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
