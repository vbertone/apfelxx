//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/qgrid.h>
#include <apfel/distribution.h>
#include <apfel/tabulateobject.h>
#include <apfel/timer.h>

#include <iostream>
#include <cmath>
#include <iomanip>

using namespace apfel;
using namespace std;

int main()
{
  // Grid
  Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,3}, SubGrid{40,8e-1,3}}, false};

  // Define distribution.
  const auto xQdist = [&](double const& x, double const& Q)->double{ return x * ( 1 - x ) * log(Q); };
  const auto d = [&](double const& Q)->Distribution{ return Distribution{g, xQdist, Q}; };

  // Tabulate distribution on a QGrid
  const TabulateObject<Distribution> dist{d, 50, 1, 1000, 3, {}};

  // Printout Qgrid
  cout << dist << endl;

  auto nx    = 10;
  auto xmin  = 1e-5;
  auto xmax  = 9e-1;
  auto xstep = exp( log( xmax / xmin ) / ( nx - 1 ) );
  auto nQ    = 5;
  auto Qmin  = 2.;
  auto Qmax  = 100.;
  auto Qstep = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  cout << "Accuracy test ..." << endl;
  cout << "Q             "
       << "x             "
       << "Interpolated  "
       << "Direct        "
       << "Ratio         " << endl;
  cout << scientific;
  auto Q = Qmin;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      auto x = xmin;
      const auto d = dist.Evaluate(Q);
      for (auto ix = 0; ix < nx; ix++)
	{
	  cout << Q << "  " << x << "  " << d.Evaluate(x) << "  " << xQdist(x,Q) << "  " << d.Evaluate(x) / xQdist(x,Q) << endl;
	  x *= xstep;
	}
      Q *= Qstep;
    }

  nx    = 1000;
  xstep = exp( log( xmax / xmin ) / ( nx - 1 ) );
  nQ    = 1000;
  Qstep = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  cout << "\nSpeed test ..." << endl;
  cout << "Interpolating " << nx << " x-space points for each of " << nQ << " Q-space points ..." << endl;
  Timer t;
  Q = Qmin;
  for (auto iQ = 0; iQ < nQ; iQ++)
    {
      auto x = xmin;
      const auto d = dist.Evaluate(Q);
      for (auto ix = 0; ix < nx; ix++)
	{
	  d.Evaluate(x);
	  x *= xstep;
	}
      Q *= Qstep;
    }
  t.stop();

  return 0;
}
