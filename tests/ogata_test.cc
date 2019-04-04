//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/ogataquadrature.h>

#include <iostream>
#include <math.h>

int main()
{
  // Ogata quadrature object
  apfel::OgataQuadrature OgataObj{};

  // Define integrand and anaylic transform
  const double a = 1;
  const auto rfunc = [=] (double const& r) -> double{ return r * exp( - pow(a * r, 2) / 2); };
  const auto kfunc = [=] (double const& k) -> double{ return exp( - pow(k / a, 2) / 2) / pow(a, 2); };

  // Print comparison
  const std::vector<double> kv = {0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6};
  std::cout << "      k         "
	    << "  analytic      "
	    << "  numerical     "
	    << "   ratio        "
	    << std::endl;
  for (auto const& k : kv)
    std::cout << std::scientific << k << "\t"
	      << kfunc(k) << "\t"
	      << OgataObj.transform(rfunc, k) << "\t"
	      << OgataObj.transform(rfunc, k) / kfunc(k) << "\t"
	      << std::endl;

  return 0;
}
