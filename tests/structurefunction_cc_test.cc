//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <cmath>
#include <map>
#include <iomanip>

int main()
{
  //SetVerbosityLevel(0);
  apfel::Banner();

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100,1e-5,3}, apfel::SubGrid{60,1e-1,3}, apfel::SubGrid{50,6e-1,3}, apfel::SubGrid{50,8e-1,3}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vectors of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 2;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // CKM matrix elements.
  std::function<std::vector<double>(double const&)> fCKM = [=] (double const&) -> std::vector<double> { return apfel::CKM2; };

  // Initialize QCD evolution objects
  const auto DglapObj = InitializeDglapObjectsQCD(g, Thresholds);

  // Construct the DGLAP object
  auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Evolved PDFs
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x,Q); };

  // Initialize coefficient functions
  const auto F2PlusCCObj  = InitializeF2CCPlusObjectsZM(g, Thresholds);
  const auto F2MinusCCObj = InitializeF2CCMinusObjectsZM(g, Thresholds);
  const auto FLPlusCCObj  = InitializeFLCCPlusObjectsZM(g, Thresholds);
  const auto FLMinusCCObj = InitializeFLCCMinusObjectsZM(g, Thresholds);
  const auto F3PlusCCObj  = InitializeF3CCPlusObjectsZM(g, Thresholds);
  const auto F3MinusCCObj = InitializeF3CCMinusObjectsZM(g, Thresholds);

  // Initialize structure functions
  const auto F2p = BuildStructureFunctions(F2PlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
  const auto F2m = BuildStructureFunctions(F2MinusCCObj, PDFs, PerturbativeOrder, as, fCKM);
  const auto FLp = BuildStructureFunctions(FLPlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
  const auto FLm = BuildStructureFunctions(FLMinusCCObj, PDFs, PerturbativeOrder, as, fCKM);
  const auto F3p = BuildStructureFunctions(F3PlusCCObj,  PDFs, PerturbativeOrder, as, fCKM);
  const auto F3m = BuildStructureFunctions(F3MinusCCObj, PDFs, PerturbativeOrder, as, fCKM);

  const apfel::TabulateObject<apfel::Distribution> F2total  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(0).Evaluate(Q) - F2m.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2light  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(1).Evaluate(Q) - F2m.at(1).Evaluate(Q) + F2p.at(2).Evaluate(Q) - F2m.at(2).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2charm  {[&] (double const& Q) -> apfel::Distribution { return F2p.at(4).Evaluate(Q) - F2m.at(4).Evaluate(Q) + F2p.at(5).Evaluate(Q) - F2m.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2bottom {[&] (double const& Q) -> apfel::Distribution { return F2p.at(3).Evaluate(Q) - F2m.at(3).Evaluate(Q) + F2p.at(6).Evaluate(Q) - F2m.at(6).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> FLtotal  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(0).Evaluate(Q) - FLm.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLlight  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(1).Evaluate(Q) - FLm.at(1).Evaluate(Q) + FLp.at(2).Evaluate(Q) - FLm.at(2).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLcharm  {[&] (double const& Q) -> apfel::Distribution { return FLp.at(4).Evaluate(Q) - FLm.at(4).Evaluate(Q) + FLp.at(5).Evaluate(Q) - FLm.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLbottom {[&] (double const& Q) -> apfel::Distribution { return FLp.at(3).Evaluate(Q) - FLm.at(3).Evaluate(Q) + FLp.at(6).Evaluate(Q) - FLm.at(6).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> F3total  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(0).Evaluate(Q) - F3m.at(0).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3light  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(1).Evaluate(Q) - F3m.at(1).Evaluate(Q) + F3p.at(2).Evaluate(Q) - F3m.at(2).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3charm  {[&] (double const& Q) -> apfel::Distribution { return F3p.at(4).Evaluate(Q) - F3m.at(4).Evaluate(Q) + F3p.at(5).Evaluate(Q) - F3m.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F3bottom {[&] (double const& Q) -> apfel::Distribution { return F3p.at(3).Evaluate(Q) - F3m.at(3).Evaluate(Q) + F3p.at(6).Evaluate(Q) - F3m.at(6).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  apfel::Timer t;
  t.start();

  // Final scale
  const double Q = 100;

  std::cout << std::scientific << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  const std::vector<double> xlha = { 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1 };

  std::cout << "    x   "
            << "  F2light   "
            << "  F2charm   "
            << "  F2bottom  "
            << "  F2total   "
            << std::endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    std::cout << std::setprecision(1) << xlha[i] << "  " << std::setprecision(4)
              << 2 * F2light.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * F2charm.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * F2bottom.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * F2total.EvaluatexQ(xlha[i],Q) << "  "
              << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
            << "  FLlight   "
            << "  FLcharm   "
            << "  FLbottom  "
            << "  FLtotal   "
            << std::endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    std::cout << std::setprecision(1) << xlha[i] << "  " << std::setprecision(4)
              << 2 * FLlight.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * FLcharm.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * FLbottom.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * FLtotal.EvaluatexQ(xlha[i],Q) << "  "
              << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
            << "  F3light   "
            << "  F3charm   "
            << "  F3bottom  "
            << "  F3total   "
            << std::endl;
  for (auto i = 2; i < (int) xlha.size(); i++)
    std::cout << std::setprecision(1) << xlha[i] << "  " << std::setprecision(4)
              << 2 * F3light.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * F3charm.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * F3bottom.EvaluatexQ(xlha[i],Q) << "  "
              << 2 * F3total.EvaluatexQ(xlha[i],Q) << "  "
              << std::endl;
  std::cout << std::endl;
  t.stop();

  const int k = 1000000;
  std::cout << "Interpolating " << k << " times F2 on the grid... ";
  t.start();
  for (auto i = 0; i < k; i++)
    F2total.EvaluatexQ(0.05,Q);
  t.stop();

  return 0;
}
