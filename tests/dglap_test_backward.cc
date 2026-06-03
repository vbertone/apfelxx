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
  const std::vector<double> Thresholds = {0, 0, 0};

  // Perturbative order
  const int PerturbativeOrder = 3;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{apfel::SubGrid{300, 1e-7, 3}, {apfel::Grid::SubGridPars{0.1, 3, 3}, apfel::Grid::SubGridPars{0.6, 3, 3}, apfel::Grid::SubGridPars{0.8, 3, 5}}};

  // DGLAP objects
  const auto DObj = InitializeDglapObjectsQCDOme(g, Thresholds);

  // Construct evolved PDF from mu0
  const auto EvolvedPDFsF = BuildDglap(DObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsF{*EvolvedPDFsF, 50, 1, 1000, 3};

  // Forward evolution
  const std::function<std::map<int, double>(double const&, double const&)> PDFsF = [=] (double const& y, double const&) -> std::map<int, double>
  {
    return TabulatedPDFsF.EvaluateMapxQ(y, mu);
  };

  // Construct evolved PDF from mu
  const auto EvolvedPDFsB = BuildDglap(DObj, PDFsF, mu, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsB{*EvolvedPDFsB, 50, 1, 1000, 3};

  // Print results
  const std::vector<double> xlha = {1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  std::cout << std::scientific;
  std::cout << "\nRatio between initial-condition PDFs at mu0 = " << mu0 << " GeV and PDFs evolved from mu0 to mu = " << mu << " GeV and back to mu0" << std::endl;
  std::cout << "   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << "  dbr-dbr   "
            << " 2(ubr+dbr) "
            << "   s+sbar   "
            << "   gluon    "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> DistMapF = apfel::QCDEvToPhys(apfel::LHToyPDFs(x, mu0));
      const std::map<int, double> DistMapB = apfel::QCDEvToPhys(TabulatedPDFsB.EvaluateMapxQ(x, mu0));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << std::setw(12) << std::right << ( DistMapF.at(2) - DistMapF.at(-2) ) / ( DistMapB.at(2) - DistMapB.at(-2) )
                << std::setw(12) << std::right << ( DistMapF.at(1) - DistMapF.at(-1) ) / ( DistMapB.at(1) - DistMapB.at(-1) )
                << std::setw(12) << std::right << ( DistMapF.at(-1) - DistMapF.at(-2) ) / ( DistMapB.at(-1) - DistMapB.at(-2) )
                << std::setw(12) << std::right << ( 2 * ( DistMapF.at(-2) + DistMapF.at(-1) ) ) / ( 2 * ( DistMapB.at(-2) + DistMapB.at(-1) ) )
                << std::setw(12) << std::right << ( DistMapF.at(3) + DistMapF.at(-3) ) / ( DistMapB.at(3) + DistMapB.at(-3) )
                << std::setw(12) << std::right << ( DistMapF.at(0) ) / ( DistMapB.at(0) )
                << std::endl;
    }

  return 0;
}
