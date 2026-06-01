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
  const int PerturbativeOrder = 1;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Get DGLAP objects in the Krk factorisation scheme
  const auto DglapObjNmKrk = ChangeFactorisationSchemeMSbarToK(InitializeDglapObjectsQCD(g, Thresholds), InitializeSchemeChangeKernelsKrk(g, Thresholds));
  const auto DglapObjAnKrk = InitializeDglapObjectsQCDKrk(g, Thresholds);

  // Construct the DGLAP objects
  const auto EvolvedPDFsNmKrk = BuildDglap(DglapObjNmKrk, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
  const auto EvolvedPDFsAnKrk = BuildDglap(DglapObjAnKrk, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsNmKrk{*EvolvedPDFsNmKrk, 50, 1, 1000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsAnKrk{*EvolvedPDFsAnKrk, 50, 1, 1000, 3};

  // Get DGLAP objects in the PHYS factorisation scheme
  const auto DglapObjNmPHYS = ChangeFactorisationSchemeMSbarToK(InitializeDglapObjectsQCD(g, Thresholds), InitializeSchemeChangeKernelsPHYS(g, Thresholds));
  const auto DglapObjAnPHYS = InitializeDglapObjectsQCDPHYS(g, Thresholds);

  // Construct the DGLAP objects
  const auto EvolvedPDFsNmPHYS = BuildDglap(DglapObjNmPHYS, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
  const auto EvolvedPDFsAnPHYS = BuildDglap(DglapObjAnPHYS, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsNmPHYS{*EvolvedPDFsNmPHYS, 50, 1, 1000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsAnPHYS{*EvolvedPDFsAnPHYS, 50, 1, 1000, 3};

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  std::cout << std::scientific;
  std::cout << "\nAlphaQCD(Q) = " << Alphas.Evaluate(mu) << std::endl;

  std::cout << "\nRatio numerical/analytic (Krk):\n   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << " 2(ubr+dbr) "
            << "   c+cbar   "
            << "    gluon   "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> DistMapNmKrk = apfel::QCDEvToPhys(TabulatedPDFsNmKrk.EvaluateMapxQ(x, mu));
      const std::map<int, double> DistMapAnKrk = apfel::QCDEvToPhys(TabulatedPDFsAnKrk.EvaluateMapxQ(x, mu));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << ( DistMapNmKrk.at(2) - DistMapNmKrk.at(-2) ) / ( DistMapAnKrk.at(2) - DistMapAnKrk.at(-2) )
                << "  " << ( DistMapNmKrk.at(1) - DistMapNmKrk.at(-1) ) / ( DistMapAnKrk.at(1) - DistMapAnKrk.at(-1) )
                << "  " << ( 2 * ( DistMapNmKrk.at(-2) + DistMapNmKrk.at(-1) ) ) / ( 2 * ( DistMapAnKrk.at(-2) + DistMapAnKrk.at(-1) ) )
                << "  " << ( DistMapNmKrk.at(4) + DistMapNmKrk.at(-4) ) / ( DistMapAnKrk.at(4) + DistMapAnKrk.at(-4) )
                << "  " << ( DistMapNmKrk.at(0) ) / ( DistMapAnKrk.at(0) )
                << std::endl;
    }

  std::cout << "\nRatio numerical/analytic (PHYS):\n   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << " 2(ubr+dbr) "
            << "   c+cbar   "
            << "    gluon   "
            << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> DistMapNmPHYS = apfel::QCDEvToPhys(TabulatedPDFsNmPHYS.EvaluateMapxQ(x, mu));
      const std::map<int, double> DistMapAnPHYS = apfel::QCDEvToPhys(TabulatedPDFsAnPHYS.EvaluateMapxQ(x, mu));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << ( DistMapNmPHYS.at(2) - DistMapNmPHYS.at(-2) ) / ( DistMapAnPHYS.at(2) - DistMapAnPHYS.at(-2) )
                << "  " << ( DistMapNmPHYS.at(1) - DistMapNmPHYS.at(-1) ) / ( DistMapAnPHYS.at(1) - DistMapAnPHYS.at(-1) )
                << "  " << ( 2 * ( DistMapNmPHYS.at(-2) + DistMapNmPHYS.at(-1) ) ) / ( 2 * ( DistMapAnPHYS.at(-2) + DistMapAnPHYS.at(-1) ) )
                << "  " << ( DistMapNmPHYS.at(4) + DistMapNmPHYS.at(-4) ) / ( DistMapAnPHYS.at(4) + DistMapAnPHYS.at(-4) )
                << "  " << ( DistMapNmPHYS.at(0) ) / ( DistMapAnPHYS.at(0) )
                << std::endl;
    }

  return 0;
}
