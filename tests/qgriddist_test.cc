//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <cmath>
#include <iomanip>

#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/qgrid.h>
#include <apfel/Distribution.h>
#include <apfel/timer.h>

using namespace apfel;
using namespace std;

/**
 * Prototype function for testing purproses.
 */
double xQdist(double const& x, double const& Q)
{
  return ( 1 - x ) * log(Q);
}

/**
 * @brief The Parton class
 */
class Function: public Distribution
{
public:
  /**
   * Allocate the langrage interpolation and fill the inherited
   * \c _distribution object with the jointed grid.
   */
  Function(Grid const& gr, double const& Q): Distribution(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(xQdist(ix,Q));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
        vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          if (ix < 1) sg.push_back(xQdist(ix,Q));
	  else        sg.push_back(0);
        _distributionSubGrid.push_back(sg);
      }
  }
};

/**
 * @brief Distibution on the grid
 */
class GridDistribution: public QGrid<Distribution>
{
public:
  GridDistribution(Grid   const& gr,
		   int    const& nQ,
		   double const& QMin,
		   double const& QMax,
		   int    const& InterDegree):
    QGrid<Distribution>(nQ, QMin, QMax, InterDegree, {})
  {
    for (auto const& iQ : _Qg)
      {
	const Function f{gr, iQ}; 
	_GridValues.push_back(f);
      }
  }

};

int main()
{
  // Grid
  Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,3}, SubGrid{40,8e-1,3}}, false};

  // Tabulate distribution on a QGrid
  const GridDistribution dist{g, 50, 1, 1000, 3};

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
  t.start();
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
  t.printTime(t.stop());

  return 0;
}
