//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <cmath>

int main()
{
  // Integrand function.
  const apfel::Integrator f{[&] (double const& x) -> double{ return log(x); }};

  // Print true value.
  std::cout << "True value: " << 2 * ( log(2) - 1 ) << std::endl;

  // Integrate using Gauss-Legendre with a given accuracy.
  const double res1 = f.integrateGL(0, 2, 1e-5);

  // Integrate using Gauss-Legendre with a different accruracy.
  const double res2 = f.integrateGL(0, 2, 1e-3);

  // Print results.
  std::cout << res1 << "  " << res2 << "  " << res1 / res2 << std::endl;

  // Integrate introducing fixed point in the integration interval.
  const double res3 = f.integrate(0, 2, {-1, 0.1, 0.5, 0.7, 1.1, 3.4}, 1e-5);

  std::cout << res1 << "  " << res3 << "  " << res1 / res3 << std::endl;

  // Now revert integration and try integrate with and without fixed
  // points.
  const double res4 = f.integrate(0, 2, 1e-5);
  const double res5 = f.integrate(0, 2, {-1, 0.1, 0.5, 0.7, 1.1, 3.4}, 1e-5);

  // Print results.
  std::cout << res4 << "  " << res5 << "  " << res4 / res5 << std::endl;

  // Integrate using Gauss-Kronrod with a given accuracy.
  const double res6 = f.integrateGK(0, 2, 1e-5);

  // Integrate using Gauss-Kronrod with a different accruracy.
  const double res7 = f.integrateGK(0, 2, 1e-3);

  // Print results.
  std::cout << res6 << "  " << res7 << "  " << res6 / res7 << std::endl;

  // Performance test
  const int k = 10000;
  std::cout << "Integrating " << k << " times with Gauss-Legendre... ";
  apfel::Timer t;
  for (int i = 0; i < k; i++)
    f.integrate(0, 2, 1e-9);
  t.stop();

  std::cout << "Integrating " << k << " times with Gauss-Kronrod... ";
  t.start();
  for (int i = 0; i < k; i++)
    f.integrate(0, 2, 1e-9);
  t.stop();

  return 0;
}
