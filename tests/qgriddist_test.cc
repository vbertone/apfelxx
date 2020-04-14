//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <cmath>
#include <iomanip>

int main()
{
  // Grid
  apfel::Grid g{{apfel::SubGrid{80,1e-5,3}, apfel::SubGrid{50,1e-1,3}, apfel::SubGrid{40,8e-1,3}}, false};

  // Thresholds
  const std::vector<double> Thrs{1.4, 6};

  // Define distribution with discontinuities at the thresholds.
  const auto xQdist = [&](double const& x, double const& Q) -> double{ return x * ( 1 - x ) / log(log(Q/0.2)) + (Q > Thrs[0] ? 0.1 : 0) + (Q > Thrs[1] ? 0.2 : 0); };

  // Define derivative in Q of the distribution.
  const auto dxQdist = [&](double const& x, double const& Q) -> double{ return - x * ( 1 - x ) / pow(log(log(Q/0.2)), 2) / Q / log(Q/0.2); };

  // Define integral in Q of the distribution.
  const auto ixQdist = [&](double const& x, double const& Qa, double const& Qb) -> double
  {
    apfel::Integrator integrand{[=] (double const& Q) -> double{ return x * ( 1 - x ) / log(log(Q/0.2)) + (Q > Thrs[0] ? 0.1 : 0) + (Q > Thrs[1] ? 0.2 : 0); }};
    return integrand.integrate(Qa, Qb, Thrs, 1e-7);
  };

  // Define distribution on the x-grid
  const auto d = [&](double const& Q) -> apfel::Distribution{ return apfel::Distribution{g, xQdist, Q}; };

  // Tabulate distribution on a QGrid
  const apfel::TabulateObject<apfel::Distribution> dist{d, 100, 1, 200, 3, Thrs};

  // Printout Qgrid
  std::cout << dist << std::endl;

  double nx    = 5;
  double xmin  = 1e-1;
  double xmax  = 9e-1;
  double xstep = exp( log( xmax / xmin ) / ( nx - 1 ) );
  double nQ    = 20;
  double Qmin  = 1.1;
  double Qmax  = 100;
  double Qstep = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  std::cout << "Accuracy test:" << std::endl;
  std::cout << "     Q        "
            << "     x        "
            << "Interpolated  "
            << "   Direct     "
            << "   Ratio      "
            << " Deriv. inter. "
            << "  Deriv. dir.  "
            << "   Ratio      "
            << "Integ. inter. "
            << " Integ. dir.  "
            << "   Ratio      "
            << std::endl;
  std::cout << std::scientific;
  double Q = Qmin;
  for (int iQ = 0; iQ < nQ; iQ++)
    {
      double x = xmin;
      const apfel::Distribution d  = dist.Evaluate(Q);
      const apfel::Distribution dd = dist.Derive(Q);
      const apfel::Distribution id = dist.Integrate(Q, 5);
      for (int ix = 0; ix < nx; ix++)
        {
          std::cout << Q << "  " << x << "  "
                    <<  d.Evaluate(x) << "  " <<  xQdist(x, Q)    << "  " <<  d.Evaluate(x) /  xQdist(x, Q) << "  "
                    << dd.Evaluate(x) << "  " << dxQdist(x, Q)    << "  " << dd.Evaluate(x) / dxQdist(x, Q) << "  "
                    << id.Evaluate(x) << "  " << ixQdist(x, Q, 5) << "  " << id.Evaluate(x) / ixQdist(x, Q, 5) << "  "
                    << std::endl;
          x *= xstep;
        }
      Q *= Qstep;
      //exit(-1);
    }

  nx    = 1000;
  xstep = exp( log( xmax / xmin ) / ( nx - 1 ) );
  nQ    = 1000;
  Qstep = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  std::cout << "\nPerformance test:" << std::endl;
  std::cout << "Interpolating " << nx << " x-space points for each of " << nQ << " Q-space points... ";
  apfel::Timer t;
  Q = Qmin;
  for (int iQ = 0; iQ < nQ; iQ++)
    {
      double x = xmin;
      const apfel::Distribution d = dist.Evaluate(Q);
      for (int ix = 0; ix < nx; ix++)
        {
          d.Evaluate(x);
          x *= xstep;
        }
      Q *= Qstep;
    }
  t.stop();

  return 0;
}
