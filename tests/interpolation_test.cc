//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <apfel/grid.h>
#include <apfel/subgrid.h>
#include <apfel/distribution.h>
#include <apfel/timer.h>

#include <iostream>
#include <iomanip>

using namespace apfel;
using namespace std;

int main()
{
  cout << setprecision(12) << scientific;

  // Grid
  const Grid g{{SubGrid{80,1e-5,3}, SubGrid{50,1e-1,5}, SubGrid{40,8e-1,5}}, false};

  // Test distribution
  const auto xg = [&] (double const& x)->double{ return x * ( 1 -x ); };
  const Distribution xgluon{g, xg};

  // Test values
  vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

  cout << "x, original function, interpolated function (joint), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix);
      cout << ix << "  "
           << original << "  "
           << interpol << "  "
           << original / interpol<< endl;
    }
  cout << endl;

  cout << "x, original function, interpolated function (first subgrid), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix,0);
      cout << ix << "  "
           << original << "  "
           << interpol << "  "
           << original / interpol<< endl;
    }
  cout << endl;

  cout << "x, original function, interpolated function (second subgrid), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix,1);
      cout << ix << "  "
           << original << "  "
           << interpol << "  "
           << original / interpol<< endl;
    }
  cout << endl;

  cout << "x, original function, interpolated function (third subgrid), ratio" << endl;
  for (auto const& ix: x)
    {
      const auto original = xg(ix);
      const auto interpol = xgluon.Evaluate(ix,2);
      cout << ix << "  "
           << original << "  "
           << interpol << "  "
           << original / interpol<< endl;
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
  t.stop();

  cout << "(First SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,0);
  t.stop();

  cout << "(Second SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,1);
  t.stop();

  cout << "(Third SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,2);
  t.stop();

  cout << "\n";

  return 0;
}
