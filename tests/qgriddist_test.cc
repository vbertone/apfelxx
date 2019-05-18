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

  // Define distribution.
  const auto xQdist = [&](double const& x, double const& Q) -> double{ return x * ( 1 - x ) * log(Q); };
  const auto d = [&](double const& Q) -> apfel::Distribution{ return apfel::Distribution{g, xQdist, Q}; };

  // Tabulate distribution on a QGrid
  const apfel::TabulateObject<apfel::Distribution> dist{d, 50, 1, 1000, 3, {}};

  // Printout Qgrid
  std::cout << dist << std::endl;

  double nx    = 10;
  double xmin  = 1e-5;
  double xmax  = 9e-1;
  double xstep = exp( log( xmax / xmin ) / ( nx - 1 ) );
  double nQ    = 5;
  double Qmin  = 2.;
  double Qmax  = 100.;
  double Qstep = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  std::cout << "Accuracy test ..." << std::endl;
  std::cout << "Q             "
            << "x             "
            << "Interpolated  "
            << "Direct        "
            << "Ratio         " << std::endl;
  std::cout << std::scientific;
  double Q = Qmin;
  for (int iQ = 0; iQ < nQ; iQ++)
    {
      double x = xmin;
      const apfel::Distribution d = dist.Evaluate(Q);
      for (int ix = 0; ix < nx; ix++)
        {
          std::cout << Q << "  " << x << "  " << d.Evaluate(x) << "  " << xQdist(x,Q) << "  " << d.Evaluate(x) / xQdist(x,Q) << std::endl;
          x *= xstep;
        }
      Q *= Qstep;
    }

  nx    = 1000;
  xstep = exp( log( xmax / xmin ) / ( nx - 1 ) );
  nQ    = 1000;
  Qstep = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  std::cout << "\nSpeed test ..." << std::endl;
  std::cout << "Interpolating " << nx << " x-space points for each of " << nQ << " Q-space points ..." << std::endl;
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
