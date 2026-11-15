//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{30, 1e-5, 3}, apfel::SubGrid{20, 1e-1, 3}, apfel::SubGrid{10, 9e-1, 3}}};
  //const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{200, 1e-1, 3}, apfel::SubGrid{100, 9e-1, 3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 1;

  // Skewness
  const double xi = 0;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Construct the DGLAP objects
  const auto EvolvedGPDs = BuildDglap(InitializeGpdObjects(g, Thresholds, xi),     apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCDOme(g, Thresholds), apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs and GPDs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 100, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*EvolvedGPDs, 30, 1, 100, 3};

  // Final scale
  const double mu = 10;

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
      const std::map<int, double> DistMapPDF = apfel::QCDEvToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu));
      const std::map<int, double> DistMapGPD = apfel::QCDEvToPhys(TabulatedGPDs.EvaluateMapxQ(x, mu));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << ( DistMapGPD.at(2) - DistMapGPD.at(-2) ) / ( DistMapPDF.at(2) - DistMapPDF.at(-2) )
                << "  " << ( DistMapGPD.at(1) - DistMapGPD.at(-1) ) / ( DistMapPDF.at(1) - DistMapPDF.at(-1) )
                << "  " << ( DistMapGPD.at(-2) + DistMapGPD.at(-1) ) / ( DistMapPDF.at(-2) + DistMapPDF.at(-1) )
                << "  " << ( DistMapGPD.at(4) + DistMapGPD.at(-4) ) / ( DistMapPDF.at(4) + DistMapPDF.at(-4) )
                << "  " << DistMapGPD.at(0) / DistMapPDF.at(0)
                << std::endl;
    }

  return 0;
}
