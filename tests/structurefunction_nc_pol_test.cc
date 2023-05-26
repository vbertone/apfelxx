//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

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

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const& Q) -> std::vector<double> { return apfel::ElectroWeakCharges(Q, false); };
  std::function<std::vector<double>(double const&)> fDq = [=] (double const& Q) -> std::vector<double> { return apfel::ParityViolatingElectroWeakCharges(Q, false); };

  // Initialize QCD evolution objects
  const auto DglapObj = InitializeDglapObjectsQCDpol(g, Thresholds);

  // Construct the DGLAP object
  auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Evolved PDFs
  const auto PDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedPDFs.EvaluateMapxQ(x, Q); };

  // Initialize coefficient functions
  const auto g4Obj = Initializeg4NCObjectsZM(g, Thresholds);
  const auto gLObj = InitializegLNCObjectsZM(g, Thresholds);
  const auto g1Obj = Initializeg1NCObjectsZM(g, Thresholds);

  // Initialize structure functions
  const auto g4 = BuildStructureFunctions(g4Obj, PDFs, PerturbativeOrder, as, fDq);
  const auto gL = BuildStructureFunctions(gLObj, PDFs, PerturbativeOrder, as, fDq);
  const auto g1 = BuildStructureFunctions(g1Obj, PDFs, PerturbativeOrder, as, fBq);

  // Tabulate Structure functions
  const apfel::TabulateObject<apfel::Distribution> g4total {[&] (double const& Q) -> apfel::Distribution{ return g4.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> g4light {[&] (double const& Q) -> apfel::Distribution{ return g4.at(1).Evaluate(Q) + g4.at(2).Evaluate(Q) + g4.at(3).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> g4charm {[&] (double const& Q) -> apfel::Distribution{ return g4.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> g4bottom{[&] (double const& Q) -> apfel::Distribution{ return g4.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> gLtotal {[&] (double const& Q) -> apfel::Distribution{ return gL.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> gLlight {[&] (double const& Q) -> apfel::Distribution{ return gL.at(1).Evaluate(Q) + gL.at(2).Evaluate(Q) + gL.at(3).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> gLcharm {[&] (double const& Q) -> apfel::Distribution{ return gL.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> gLbottom{[&] (double const& Q) -> apfel::Distribution{ return gL.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

  const apfel::TabulateObject<apfel::Distribution> g1total {[&] (double const& Q) -> apfel::Distribution{ return g1.at(0).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> g1light {[&] (double const& Q) -> apfel::Distribution{ return g1.at(1).Evaluate(Q) + g1.at(2).Evaluate(Q) + g1.at(3).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> g1charm {[&] (double const& Q) -> apfel::Distribution{ return g1.at(4).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> g1bottom{[&] (double const& Q) -> apfel::Distribution{ return g1.at(5).Evaluate(Q); }, 50, 1, 200, 3, Thresholds};

  apfel::Timer t;

  // Final scale
  const double Q = 100;

  std::cout << std::scientific << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << "    x   "
            << "  g4light   "
            << "  g4charm   "
            << "  g4bottom  "
            << "  g4total   "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(1) << x << "  " << std::setprecision(4)
              << g4light.EvaluatexQ(x, Q)  << "  "
              << g4charm.EvaluatexQ(x, Q)  << "  "
              << g4bottom.EvaluatexQ(x, Q) << "  "
              << g4total.EvaluatexQ(x, Q)  << "  "
              << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
            << "  gLlight   "
            << "  gLcharm   "
            << "  gLbottom  "
            << "  gLtotal   "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(1) << x << "  " << std::setprecision(4)
              << gLlight.EvaluatexQ(x, Q)  << "  "
              << gLcharm.EvaluatexQ(x, Q)  << "  "
              << gLbottom.EvaluatexQ(x, Q) << "  "
              << gLtotal.EvaluatexQ(x, Q)  << "  "
              << std::endl;
  std::cout << std::endl;

  std::cout << "    x   "
            << "  g1light   "
            << "  g1charm   "
            << "  g1bottom  "
            << "  g1total   "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(1) << x << "  " << std::setprecision(4)
              << g1light.EvaluatexQ(x, Q)  << "  "
              << g1charm.EvaluatexQ(x, Q)  << "  "
              << g1bottom.EvaluatexQ(x, Q) << "  "
              << g1total.EvaluatexQ(x, Q)  << "  "
              << std::endl;
  std::cout << std::endl;

  t.stop();

  const int k = 1000000;
  std::cout << "Interpolating " << k << " times g1 on the grid... ";
  t.start();
  for (int i = 0; i < k; i++)
    g1total.EvaluatexQ(0.05, Q);
  t.stop();

  return 0;
}
