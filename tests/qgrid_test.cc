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
  // Constructor of QGrid with type double
  const apfel::QGrid<double> qg{50, 1, 1000, 3, {0, 0, 0, sqrt(2), 4.5, 175.}};
  std::cout << qg << std::endl;

  // AlphaQCD running coupling
  apfel::AlphaQCD as{0.35, sqrt(2), {0, 0, 0, sqrt(2), 4.5, 175}, 2};

  // Tabulate AlphaQCD on a QGrid
  const apfel::TabulateObject<double> gas{as, 50, 1, 1000, 3};

  std::cout << "\nAccuracy test:" << std::endl;
  double nQ   = 20;
  double Qmin = 1.1;
  double Qmax = 999.;
  double Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );
  double Q = Qmin;
  std::cout << std::setprecision(8) << std::scientific;
  std::cout << "Q       \t\tDirect  \t\tInterpolated\t\tRatio" << std::endl;
  for (int iQ = 0; iQ < nQ; iQ++)
    {
      std::cout << Q << "\t\t" << as.Evaluate(Q) << "\t\t" << gas.Evaluate(Q) << "\t\t" << as.Evaluate(Q) / gas.Evaluate(Q) << std::endl;
      Q *=Step;
    }

  std::cout << "\nPerformance test:" << std::endl;
  apfel::Timer t;
  nQ   = 100000;
  Step = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  t.start();
  std::cout << "Direct calculation of " << nQ << " points... ";
  Q = Qmin;
  for (int iQ = 0; iQ < nQ; iQ++)
    {
      as.Evaluate(Q);
      Q *= Step;
    }
  t.stop();

  t.start();
  std::cout << "Interpolated calculation of " << nQ << " points... ";
  Q = Qmin;
  for (int iQ = 0; iQ < nQ; iQ++)
    {
      gas.Evaluate(Q);
      Q *= Step;
    }
  t.stop();

  // Test behavior of qgrid for Q in the gap between thresholds
  const std::vector<double> thresholds{0, 0, 0, 2., 4.5, 175.};
  const auto count = [&thresholds](double const& Q) -> double { return apfel::NF(Q, thresholds); };
  const apfel::TabulateObject<double> table{count, 100, 1., 10., 3, thresholds};
  const auto& grid = table.GetQGrid();
  for (const double threshold : {2., 4.5})
    {
      const auto upper = std::upper_bound(grid.begin(), grid.end(), threshold);
      const double below = *(upper - 1);
      const double above = *upper;
      for (const double Q : {below, std::nextafter(threshold, 0.), threshold, std::nextafter(threshold, std::numeric_limits<double>::infinity()),above})
	if (std::abs(table.Evaluate(Q) - count(Q)) > 1.e-12)
	  throw std::runtime_error("QGrid threshold-gap evaluation disagrees with NF");
    }

    return 0;
}
