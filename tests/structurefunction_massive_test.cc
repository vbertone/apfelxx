//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
  const std::vector<double> Thresholds{0, 0, 0};
  const std::vector<double> Masses{0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 1;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };

  // Initialize QCD evolution objects
  const auto DglapObj = InitializeDglapObjectsQCD(g, Thresholds);

  // Construct the DGLAP object
  auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Evolved PDFs
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int,double> { return TabulatedPDFs.EvaluateMapxQ(x,Q); };

  // Initialize coefficient functions
  const auto F2Obj = InitializeF2NCObjectsMassive(g, Masses);
  const auto FLObj = InitializeFLNCObjectsMassive(g, Masses);

  // Initialize structure functions
  const auto F2 = BuildStructureFunctions(F2Obj, PDFs, PerturbativeOrder, as, fBq);
  const auto FL = BuildStructureFunctions(FLObj, PDFs, PerturbativeOrder, as, fBq);

  // Tabulate Structure functions
  const apfel::TabulateObject<apfel::Distribution> F2charm {[&] (double const& Q) -> apfel::Distribution{ return F2.at(4).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> F2bottom{[&] (double const& Q) -> apfel::Distribution{ return F2.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLcharm {[&] (double const& Q) -> apfel::Distribution{ return FL.at(4).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> FLbottom{[&] (double const& Q) -> apfel::Distribution{ return FL.at(5).Evaluate(Q); }, 50, 1, 1000, 3, Thresholds};

  apfel::Timer t;

  // Final scale
  const double Q = 5;

  // Scaling variables for charm and bottom
  const double etac = Q * Q / ( Q * Q + 4 * Masses[3] * Masses[3] );
  const double etab = Q * Q / ( Q * Q + 4 * Masses[4] * Masses[4] );

  // Print results
  std::cout << std::scientific << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  std::cout << "    x   "
            << "  F2charm   "
            << "  F2bottom  "
            << std::endl;
  for (auto const& x : xlha)
    std::cout << std::setprecision(1) << x << "  " << std::setprecision(4)
              << F2charm.EvaluatexQ(x / etac, Q)  << "  "
              << F2bottom.EvaluatexQ(x / etab, Q)  << "  "
              << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
            << "  FLcharm   "
            << "  FLbottom  "
            << std::endl;
  for (auto const& x : xlha)
    std::cout << std::setprecision(1) << x << "  " << std::setprecision(4)
              << FLcharm.EvaluatexQ(x / etac, Q)  << "  "
              << FLbottom.EvaluatexQ(x / etab, Q)  << "  "
              << std::endl;
  std::cout << std::endl;
  t.stop();

  return 0;
}
