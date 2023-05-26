//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
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
  const apfel::Grid g{{apfel::SubGrid{200, 1e-9, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Get PDFs at the final scale as distributions
  const apfel::Set<apfel::Distribution> Dists = TabulatedPDFs.Evaluate(mu);

  // Momentum sum rule
  std::cout << "\nMomentum sum rule: " << Dists.at(0).Integrate(1e-9, 1) + Dists.at(1).Integrate(1e-9, 1) << std::endl;

  // Total valence sum rule
  std::cout << "Total valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * Dists.at(2)).Integrate(1e-9, 1) << std::endl;

  // V3 valence sum rule
  std::cout << "V3 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * Dists.at(4)).Integrate(1e-9, 1) << std::endl;

  // V8 valence sum rule
  std::cout << "V8 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * Dists.at(6)).Integrate(1e-9, 1) << std::endl;

  // Rotate PDFs into the physical basis
  const std::map<int, apfel::Distribution> DistsPhys = apfel::QCDEvToPhys(Dists.GetObjects());

  // Up valence sum rule
  std::cout << "Up valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(2) - DistsPhys.at(-2) )).Integrate(1e-9, 1) << std::endl;

  // Down valence sum rule
  std::cout << "Down valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(1) - DistsPhys.at(-1) )).Integrate(1e-9, 1) << std::endl;

  // Strange valence sum rule
  std::cout << "Strange valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(3) - DistsPhys.at(-3) )).Integrate(1e-9, 1) << std::endl;

  return 0;
}
