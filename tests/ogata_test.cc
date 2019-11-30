//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <math.h>

int main()
{
  // Ogata quadrature object
  const apfel::OgataQuadrature OgataObj{};

  // Double exponential quadrature object
  const apfel::DoubleExponentialQuadrature DEObj{1e-7};

  // Define integrand and anaylic transform
  const double a = 1;
  const std::function<double(double const&)> rfunc = [=] (double const& r) -> double{ return r * exp( - pow(a * r, 2) / 2); };
  const std::function<double(double const&)> kfunc = [=] (double const& k) -> double{ return exp( - pow(k / a, 2) / 2) / pow(a, 2); };

  // Print comparison
  const std::vector<double> kv = {0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6};
  std::cout << "     k          "
            << " analytic       "
            << "numerical Og    "
            << "numerical DE    "
            << "  ratio Og      "
            << "  ratio DE      "
            << std::endl;
  for (auto const& k : kv)
    std::cout << std::scientific << k << "\t"
              << kfunc(k) << "\t"
              << OgataObj.transform(rfunc, k) << "\t"
              << DEObj.transform(rfunc, k) << "\t"
              << OgataObj.transform(rfunc, k) / kfunc(k) << "\t"
              << DEObj.transform(rfunc, k) / kfunc(k) << "\t"
              << std::endl;

  return 0;
}
