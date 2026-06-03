//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iomanip>

int main()
{
  //SetVerbosityLevel(apfel::verbosity::LOW);
  apfel::Banner();

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vectors of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 1;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };
  std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return apfel::ParityViolatingElectroWeakCharges(Q, false); };

  // Initialize QCD evolution objects
  const auto KKrk = InitializeSchemeChangeKernelsKrk(g, Thresholds);
  const auto DglapObjMSb = InitializeDglapObjectsQCD(g, Thresholds);
  const auto DglapObjKrk = ChangeFactorisationSchemeMSbarToK(DglapObjMSb, KKrk);

  // Construct the DGLAP object
  auto EvolvedPDFsMSb = BuildDglap(DglapObjMSb, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);
  auto EvolvedPDFsKrk = BuildDglap(DglapObjKrk, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsMSb{*EvolvedPDFsMSb, 50, 1, 1000, 3};
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFsKrk{*EvolvedPDFsKrk, 50, 1, 1000, 3};

  // Evolved PDFs
  const auto PDFsMSb = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFsMSb.EvaluateMapxQ(x, Q); };
  const auto PDFsKrk = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFsKrk.EvaluateMapxQ(x, Q); };

  // Initialize coefficient functions
  const auto F2ObjMSb = InitializeF2NCObjectsZM(g, Thresholds);
  const auto FLObjMSb = InitializeFLNCObjectsZM(g, Thresholds);
  const auto F3ObjMSb = InitializeF3NCObjectsZM(g, Thresholds);
  const auto F2ObjKrk = ChangeFactorisationSchemeMSbarToK(F2ObjMSb, KKrk);
  const auto FLObjKrk = ChangeFactorisationSchemeMSbarToK(FLObjMSb, KKrk);
  const auto F3ObjKrk = ChangeFactorisationSchemeMSbarToK(F3ObjMSb, KKrk);

  // Initialize structure functions
  const auto F2MSb = BuildStructureFunctions(F2ObjMSb, PDFsMSb, PerturbativeOrder, as, fBq);
  const auto FLMSb = BuildStructureFunctions(FLObjMSb, PDFsMSb, PerturbativeOrder, as, fBq);
  const auto F3MSb = BuildStructureFunctions(F3ObjMSb, PDFsMSb, PerturbativeOrder, as, fDq);

  const auto F2Krk = BuildStructureFunctions(F2ObjKrk, PDFsKrk, PerturbativeOrder, as, fBq);
  const auto FLKrk = BuildStructureFunctions(FLObjKrk, PDFsKrk, PerturbativeOrder, as, fBq);
  const auto F3Krk = BuildStructureFunctions(F3ObjKrk, PDFsKrk, PerturbativeOrder, as, fDq);

  // Tabulate Structure functions
  const apfel::TabulateObject<apfel::Distribution> F2totalMSb {[&] (double const& Q) -> apfel::Distribution{ return F2MSb.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotalMSb {[&] (double const& Q) -> apfel::Distribution{ return FLMSb.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3totalMSb {[&] (double const& Q) -> apfel::Distribution{ return F3MSb.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> F2totalKrk {[&] (double const& Q) -> apfel::Distribution{ return F2Krk.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLtotalKrk {[&] (double const& Q) -> apfel::Distribution{ return FLKrk.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3totalKrk {[&] (double const& Q) -> apfel::Distribution{ return F3Krk.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

  apfel::Timer t;

  // Final scale
  const double Q = 100;

  std::cout << std::scientific << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << "    x   "
            << "  F2(MSb)   "
            << "  F2(Krk)   "
            << "  FL(MSb)   "
            << "  FL(Krk)   "
            << "  F3(MSb)   "
            << "  F3(Krk)   "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(1) << x << "  " << std::setprecision(4)
              << F2totalMSb.EvaluatexQ(x, Q)  << "  "
              << F2totalKrk.EvaluatexQ(x, Q)  << "  "
              << FLtotalMSb.EvaluatexQ(x, Q)  << "  "
              << FLtotalKrk.EvaluatexQ(x, Q)  << "  "
              << F3totalMSb.EvaluatexQ(x, Q)  << "  "
              << F3totalKrk.EvaluatexQ(x, Q)  << "  "
              << std::endl;
  t.stop();
}
