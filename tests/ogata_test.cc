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
  const apfel::DoubleExponentialQuadrature DEObj{};

  // Define integrand and anaylic transform
  const double a = 1;
  const std::function<double(double const&)> rfunc = [=] (double const& r) -> double{ return r * exp( - pow(a * r, 2) / 2); };
  const std::function<double(double const&)> kfunc = [=] (double const& k) -> double{ return exp( - pow(k / a, 2) / 2) / pow(a, 2); };

  // Print comparison
  const std::vector<double> kv = {0.0001, 0.001, 0.1, 0.2, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6};
  std::cout << "\nAccuracy test:" << std::endl;
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

  std::cout << "\nPerformance test:" << std::endl;
  apfel::Timer t;
  const int nqT = 1000000;

  t.start();
  std::cout << "Ogata quadrature: calculation of " << nqT << " trasforms... ";
  for (int iqT = 0; iqT < nqT; iqT++)
    OgataObj.transform(rfunc, 1);
  t.stop();

  t.start();
  std::cout << "Double exponential quadrature: calculation of " << nqT << " trasforms... ";
  for (int iqT = 0; iqT < nqT; iqT++)
    DEObj.transform(rfunc, 1);
  t.stop();

  // Integrands for the quadrature
  const std::function<double(double const&)> rfuncI = [=] (double const& r) -> double{ return rfunc(r) / r; };
  const std::function<double(double const&)> kfuncI = [=] (double const& k) -> double{ return k * kfunc(k); };

  // Define integrand
  const apfel::Integrator NumInt{kfuncI};

  // Ogata quadrature object for the integral
  const apfel::OgataQuadrature OgataObjI{1, apfel::eps7};

  // Double exponential quadrature object for the integral
  const apfel::DoubleExponentialQuadrature DEObjI{1};

  std::cout << "\nIntegration test:" << std::endl;
  std::cout << "   kmin         "
            << "   kmax         "
            << " Numerical      "
            << "Analytic Og     "
            << "Analytic DE     "
            << "  ratio Og      "
            << "  ratio DE      "
            << std::endl;
  for (int ik = 1; ik < (int) kv.size() - 1; ik++)
    {
      const double INum = NumInt.integrate(kv[ik], kv[ik+1], apfel::eps5);
      const double IOga = kv[ik+1] * OgataObjI.transform(rfuncI, kv[ik+1]) - kv[ik] * OgataObjI.transform(rfuncI, kv[ik]);
      const double IDEQ = kv[ik+1] * DEObjI.transform(rfuncI, kv[ik+1]) - kv[ik] * DEObjI.transform(rfuncI, kv[ik]);

      std::cout << std::scientific << kv[ik] << "\t" << kv[ik+1] << "\t"
                << INum << "\t"
                << IOga << "\t"
                << IDEQ << "\t"
                << IOga / INum << "\t"
                << IDEQ / INum << "\t"
                << std::endl;
    }

  return 0;
}
