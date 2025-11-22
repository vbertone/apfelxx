//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include <iomanip>

int main()
{
  apfel::Banner();

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{200, 1e-11, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 0;

  // Skewness
  const double xi = 0.27;

  // Running coupling
  const double AlphaQCDRef = 0.35;
  const double MuAlphaQCDRef = sqrt(2);
  apfel::AlphaQCD a{AlphaQCDRef, MuAlphaQCDRef, Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // Effective charges
  std::function<std::vector<double>(double const&)> fBq = [=] (double const&) -> std::vector<double> { return apfel::QCh2; };

  // Initialize GPD evolution objects
  const auto GpdObj = InitializeGpdObjects(g, Thresholds, xi);

  // Define set of toy GPDs
  const auto ToyGPDs = [=] (double const& x, double const&) -> std::map<int, double>
  {
    const double gl  = x * ( 1 - x * x ) / x;
    const double up  = x * x * ( 1 - x * x );
    const double upb = - x * x * x * ( 1 - x * x );
    const double dn  = x * ( 1 - x * x ) * std::exp(- pow(x - xi, 2));
    const double dnb = - x * ( 1 - x * x ) * std::exp(- pow(- x - xi, 2));

    // Construct QCD evolution basis conbinations
    double const Gluon   = gl;
    double const Singlet = dn + dnb + up + upb;
    double const T3      = up + upb - dn - dnb;
    double const Valence = dn - dnb + up - upb;
    double const V3      = up - upb - dn + dnb;

    // Fill in map in the QCD evolution basis
    std::map<int, double> QCDEvMap;
    QCDEvMap[0]  = Gluon;
    QCDEvMap[1]  = Singlet;
    QCDEvMap[2]  = Valence;
    QCDEvMap[3]  = T3;
    QCDEvMap[4]  = V3;
    QCDEvMap[5]  = Singlet;
    QCDEvMap[6]  = Valence;
    QCDEvMap[7]  = Singlet;
    QCDEvMap[8]  = Valence;
    QCDEvMap[9]  = Singlet;
    QCDEvMap[10] = Valence;
    QCDEvMap[11] = Singlet;
    QCDEvMap[12] = Valence;

    return QCDEvMap;
  };

  // Construct the DGLAP objects
  const auto EvolvedGPDs = BuildDglap(GpdObj, ToyGPDs, mu0, PerturbativeOrder, as);

  // Tabulate GPDs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedGPDs{*EvolvedGPDs, 30, 1, 100, 3};

  // Evolved GPDs
  const auto GPDs = [&] (double const& x, double const& Q) -> std::map<int, double> { return TabulatedGPDs.EvaluateMapxQ(x, Q); };

  // Initialize coefficient functions for CFFs
  const auto ImF1Obj = InitializeImCFF1NCObjectsZM(g, Thresholds, xi);
  const auto ReF1Obj = InitializeReCFF1NCObjectsZM(g, Thresholds, xi);

  // Initialize form factors
  const auto ImF1 = BuildStructureFunctions(ImF1Obj, GPDs, PerturbativeOrder, as, fBq);
  const auto ReF1 = BuildStructureFunctions(ReF1Obj, GPDs, PerturbativeOrder, as, fBq);

  // Tabulate Structure functions
  const apfel::TabulateObject<apfel::Distribution> ImF1total {[&] (double const& Q) -> apfel::Distribution{ return ImF1.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};
  const apfel::TabulateObject<apfel::Distribution> ReF1total {[&] (double const& Q) -> apfel::Distribution{ return ReF1.at(0).Evaluate(Q); }, 50, 1, 99, 3, Thresholds};

  apfel::Timer t;

  // Final scale
  const double Q = 10;

  std::cout << std::scientific << std::endl;
  std::cout << "xi = " << xi << std::endl;
  std::cout << "Alphas(Q) = " << as(Q) << std::endl;
  std::cout << std::endl;

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};

  std::cout << "     x    "
            << "     ReF1tot    "
            << "     ImF1tot    "
            << std::endl;
  for (double const& x : xlha)
    std::cout << std::setprecision(4) << x << "  " << std::setprecision(8)
              << ReF1total.EvaluatexQ(x, Q) / x << "  "
              << ImF1total.EvaluatexQ(x, Q) / x << "  "
              << std::endl;
  std::cout << std::endl;

  t.stop();

  const int k = 1000000;
  std::cout << "Interpolating " << k << " times F1 on the grid... ";
  t.start();
  for (int i = 0; i < k; i++)
    ReF1total.EvaluatexQ(0.05, Q);
  t.stop();

  return 0;
}
