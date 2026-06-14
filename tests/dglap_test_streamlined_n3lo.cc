//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include<iomanip>

int main()
{
  apfel::Banner();

  // Initial scale
  const double mu0 = sqrt(2);

  // Final scale
  const double mu = 100;

  // Vectors of masses and thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 3;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{apfel::SubGrid{300, 1e-7, 3}, {apfel::Grid::SubGridPars{0.1, 3, 3}, apfel::Grid::SubGridPars{0.6, 3, 3}, apfel::Grid::SubGridPars{0.8, 3, 5}}};

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCDOme(g, Thresholds), apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Print results
  const std::vector<double> xlha = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  std::cout << std::scientific;
  std::cout << "\nAlphaQCD(Q) = " << Alphas.Evaluate(mu) << std::endl;
  std::cout << "\n   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << "  dbr-ubr   "
            << " 2(ubr+dbr) "
            << "   s-sbar   "
            << "   s+sbar   "
            << "   c+cbar   "
            << "   b+bbar   "
            << "   gluon    "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> DistMap = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << std::setw(12) << std::right << DistMap.at(2) - DistMap.at(-2)
                << std::setw(12) << std::right << DistMap.at(1) - DistMap.at(-1)
                << std::setw(12) << std::right << DistMap.at(-1) - DistMap.at(-2)
                << std::setw(12) << std::right << 2 * ( DistMap.at(-2) + DistMap.at(-1) )
                << std::setw(12) << std::right << DistMap.at(3) - DistMap.at(-3)
                << std::setw(12) << std::right << DistMap.at(3) + DistMap.at(-3)
                << std::setw(12) << std::right << DistMap.at(4) + DistMap.at(-4)
                << std::setw(12) << std::right << DistMap.at(5) + DistMap.at(-5)
                << std::setw(12) << std::right << DistMap.at(0)
                << std::endl;
    }

  return 0;
}
