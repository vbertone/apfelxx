//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iostream>
#include <iomanip>

int main()
{
  std::cout << std::setprecision(12) << std::scientific;

  // Grid
  const apfel::Grid g{{apfel::SubGrid{80,1e-5,3}, apfel::SubGrid{50,1e-1,5}, apfel::SubGrid{40,8e-1,5}}, false};

  // Test distribution
  const auto xg = [&] (double const& x) -> double{ return x * ( 1 - x ); };
  const apfel::Distribution xgluon{g, xg};

  // Test values
  std::vector<double> x = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9};

  std::cout << "x, original function, interpolated function (joint), ratio" << std::endl;
  for (auto const& ix: x)
    {
      const double original = xg(ix);
      const double interpol = xgluon.Evaluate(ix);
      std::cout << ix << "  "
		<< original << "  "
		<< interpol << "  "
		<< original / interpol<< std::endl;
    }
  std::cout << "\n";

  std::cout << "x, original function, interpolated function (first subgrid), ratio" << std::endl;
  for (auto const& ix: x)
    {
      const double original = xg(ix);
      const double interpol = xgluon.Evaluate(ix,0);
      std::cout << ix << "  "
		<< original << "  "
		<< interpol << "  "
		<< original / interpol<< std::endl;
    }
  std::cout << "\n";

  std::cout << "x, original function, interpolated function (second subgrid), ratio" << std::endl;
  for (auto const& ix: x)
    {
      const double original = xg(ix);
      const double interpol = xgluon.Evaluate(ix,1);
      std::cout << ix << "  "
		<< original << "  "
		<< interpol << "  "
		<< original / interpol<< std::endl;
    }
  std::cout << "\n";

  std::cout << "x, original function, interpolated function (third subgrid), ratio" << std::endl;
  for (auto const& ix: x)
    {
      const double original = xg(ix);
      const double interpol = xgluon.Evaluate(ix,2);
      std::cout << ix << "  "
		<< original << "  "
		<< interpol << "  "
		<< original / interpol<< std::endl;
    }
  std::cout << "\n";

  const int nint = 1000000;
  const apfel::SubGrid test_grid{nint, 1e-5, 1};
  apfel::Timer t;

  std::cout << "Performance test ("<< nint << " interpolations) ..." << std::endl;

  std::cout << "(Joint Grid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r);
  t.stop();

  std::cout << "(First SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,0);
  t.stop();

  std::cout << "(Second SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,1);
  t.stop();

  std::cout << "(Third SubGrid) ";
  t.start();
  for (auto const& r: test_grid.GetGrid())
    xgluon.Evaluate(r,2);
  t.stop();

  std::cout << "\n";

  return 0;
}
