//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

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
  const int PerturbativeOrder = 2;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{100,5e-1,3}, apfel::SubGrid{80,8e-1,5}}};

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  std::cout << std::scientific;
  std::cout << "\nAlphaQCD(Q) = " << Alphas.Evaluate(mu) << std::endl;
  std::cout << "\n   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << " 2(ubr+dbr) "
            << "   c+cbar   "
            << "    gluon   "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> DistMap = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x,mu));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << DistMap.at(2) - DistMap.at(-2)
                << "  " << DistMap.at(1) - DistMap.at(-1)
                << "  " << 2 * ( DistMap.at(-2) + DistMap.at(-1) )
                << "  " << DistMap.at(4) + DistMap.at(-4)
                << "  " << DistMap.at(0)
                << std::endl;
    }

  return 0;
}
