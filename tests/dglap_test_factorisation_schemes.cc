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
  const double mu = 999999;

  // Vectors of masses and thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 1;

  // Lower integration limit
  const double xmin = 1e-9;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 200, 0.9, 1000001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{200, 1e-9, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Sets of input PDFs in the MSbar scheme
  const apfel::Set<apfel::Distribution> SetLHToyPDFs{apfel::EvolutionBasisQCD{apfel::NF(mu0, Thresholds)}, DistributionMap(g, apfel::LHToyPDFs, mu0)};

  // Get DGLAP objects in the Krk factorisation scheme
  const auto KKrk = InitializeSchemeChangeKernelsKrk(g, Thresholds);
  const auto DglapObjNmKrk = ChangeFactorisationSchemeMSbarToK(InitializeDglapObjectsQCD(g, Thresholds), KKrk);
  const auto DglapObjAnKrk = InitializeDglapObjectsQCDKrk(g, Thresholds);

  // Trasformed sets of input PDFs
  const apfel::Set<apfel::Distribution> SetLHToyPDFsKrk = ChangeFactorisationSchemeMSbarToK(SetLHToyPDFs, KKrk, apfel::NF(mu0, Thresholds), as(mu0));
  const std::function<std::map<int, double>(double const&, double const&)> LHToyPDFsKrk = [=] (double const& x, double const&) -> std::map<int, double>
  {
    std::map<int, double> out;
    for (auto const& d : SetLHToyPDFsKrk.GetObjects())
      out.insert({d.first, d.second.Evaluate(x)});
    return out;
  };

  // Construct the DGLAP objects
  const auto EvolvedPDFsNmKrk = BuildDglap(DglapObjNmKrk, LHToyPDFsKrk, mu0, PerturbativeOrder, as);
  const auto EvolvedPDFsAnKrk = BuildDglap(DglapObjAnKrk, LHToyPDFsKrk, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsNmKrk{*EvolvedPDFsNmKrk, 200, 1, 1000000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsAnKrk{*EvolvedPDFsAnKrk, 200, 1, 1000000, 3};

  // Get DGLAP objects in the PHYS factorisation scheme
  const auto KPHYS = InitializeSchemeChangeKernelsPHYS(g, Thresholds);
  const auto DglapObjNmPHYS = ChangeFactorisationSchemeMSbarToK(InitializeDglapObjectsQCD(g, Thresholds), KPHYS);
  const auto DglapObjAnPHYS = InitializeDglapObjectsQCDPHYS(g, Thresholds);

  // Trasformed sets of input PDFs
  // Trasformed sets of input PDFs
  const apfel::Set<apfel::Distribution> SetLHToyPDFsPHYS = ChangeFactorisationSchemeMSbarToK(SetLHToyPDFs, KPHYS, apfel::NF(mu0, Thresholds), as(mu0));
  const std::function<std::map<int, double>(double const&, double const&)> LHToyPDFsPHYS = [=] (double const& x, double const&) -> std::map<int, double>
  {
    std::map<int, double> out;
    for (auto const& d : SetLHToyPDFsPHYS.GetObjects())
      out.insert({d.first, d.second.Evaluate(x)});
    return out;
  };

  // Construct the DGLAP objects
  const auto EvolvedPDFsNmPHYS = BuildDglap(DglapObjNmPHYS, LHToyPDFsPHYS, mu0, PerturbativeOrder, as);
  const auto EvolvedPDFsAnPHYS = BuildDglap(DglapObjAnPHYS, LHToyPDFsPHYS, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsNmPHYS{*EvolvedPDFsNmPHYS, 200, 1, 1000000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsAnPHYS{*EvolvedPDFsAnPHYS, 200, 1, 1000000, 3};

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

  // Sum rules in the Krk scheme
  const apfel::Set<apfel::Distribution> DistsKrk = TabulatedPDFsNmKrk.Evaluate(mu);
  std::cout << "\nMomentum sum rule: " << DistsKrk.at(0).Integrate(xmin, 1) + DistsKrk.at(1).Integrate(xmin, 1) << std::endl;
  std::cout << "Total valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * DistsKrk.at(2)).Integrate(xmin, 1) << std::endl;
  std::cout << "V3 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * DistsKrk.at(4)).Integrate(xmin, 1) << std::endl;
  std::cout << "V8 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * DistsKrk.at(6)).Integrate(xmin, 1) << std::endl;
  const std::map<int, apfel::Distribution> DistsKrkP = QCDEvToPhys(DistsKrk.GetObjects());
  std::cout << "d sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsKrkP.at(1) - DistsKrkP.at(-1) )).Integrate(xmin, 1) << std::endl;
  std::cout << "u sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsKrkP.at(2) - DistsKrkP.at(-2) )).Integrate(xmin, 1) << std::endl;

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

  // Sum rules in the PHYS scheme
  const apfel::Set<apfel::Distribution> DistsPHYS = TabulatedPDFsNmPHYS.Evaluate(mu);
  std::cout << "\nMomentum sum rule: " << DistsPHYS.at(0).Integrate(xmin, 1) + DistsPHYS.at(1).Integrate(xmin, 1) << std::endl;
  std::cout << "Total valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * DistsPHYS.at(2)).Integrate(xmin, 1) << std::endl;
  std::cout << "V3 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * DistsPHYS.at(4)).Integrate(xmin, 1) << std::endl;
  std::cout << "V8 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * DistsPHYS.at(6)).Integrate(xmin, 1) << std::endl;
  const std::map<int, apfel::Distribution> DistsPHYSP = QCDEvToPhys(DistsPHYS.GetObjects());
  std::cout << "d sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPHYSP.at(1) - DistsPHYSP.at(-1) )).Integrate(xmin, 1) << std::endl;
  std::cout << "u sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPHYSP.at(2) - DistsPHYSP.at(-2) )).Integrate(xmin, 1) << std::endl;

  return 0;
}
