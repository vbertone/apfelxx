//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <iomanip>
#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/interpolator.h>
#include <apfel/lagrangeinterpolator.h>
#include <apfel/timer.h>
using namespace apfel;
using namespace std;

/**
 * Prototype function for testing purproses.
 */
double xg(double const& x)
{
  return x * ( 1 - x );
}

/**
 * @brief The Parton class
 */
class Parton: public LagrangeInterpolator
{
public:
  /**
   * Allocate the langrage interpolation and fill the inherited
   * \c _distribution object with the jointed grid.
   */
  Parton(Grid const& gr): LagrangeInterpolator(gr)
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      if (ix < 1) _distributionJointGrid.push_back(xg(ix));
      else        _distributionJointGrid.push_back(0);

    for (auto ig=0; ig<_grid.nGrids(); ig++)
      {
	vector<double> sg;
	for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
	  if (ix < 1) sg.push_back(xg(ix));
	  else        sg.push_back(0);
	_distributionSubGrid.push_back(sg);
      }
  }
};

int main()
{
  cout << setprecision(15) << scientific;

  const Grid g{
    {SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false
  };

  const Parton xgluon{g};

  vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

  cout << "x, original function, interpolated function (joint), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix);
      cout << ix << " "
           << original << " "
           << interpol << " "
           << original/interpol<< endl;
    }
  cout << endl;

  cout << "x, original function, interpolated function (first subgrid), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix,0);
      cout << ix << " "
           << original << " "
           << interpol << " "
           << original/interpol<< endl;
    }
  cout << endl;

  cout << "x, original function, interpolated function (second subgrid), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix,1);
      cout << ix << " "
           << original << " "
           << interpol << " "
           << original/interpol<< endl;
    }
  cout << endl;

  cout << "x, original function, interpolated function (third subgrid), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix,2);
      cout << ix << " "
           << original << " "
           << interpol << " "
           << original/interpol<< endl;
    }
  cout << endl;

  int const nint = 1000000;
  const SubGrid test_grid{nint, 1e-5, 1};
  Timer t;

  cout << "Performance test ("<< nint << " interpolations) ..." << endl;

  cout << "(Joint Grid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r);
  t.printTime(t.stop());

  cout << "(First SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,0);
  t.printTime(t.stop());

  cout << "(Second SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,1);
  t.printTime(t.stop());

  cout << "(Third SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,2);
  t.printTime(t.stop());

  cout << "\n";

  return 0;
}
