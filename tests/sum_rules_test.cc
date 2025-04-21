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
  const int PerturbativeOrder = 3;

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  const apfel::TabulateObject<double> Alphas{a, 100, 0.9, 1001, 3};
  const auto as = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{200, 1e-9, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Initialise DGLAP objects
  const auto DglapObj = InitializeDglapObjectsQCD(g, Thresholds, false, apfel::eps5, true, {0, 0, 0, 0, 0, 0, 0});

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFs, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Get PDFs at the final scale as distributions
  const apfel::Set<apfel::Distribution> Dists = TabulatedPDFs.Evaluate(mu);

  // Lower integration limit
  const double xmin = 1e-9;

  // Start the checks
  std::cout << "\nChecking sum rules at the level of PDFs..." << std::endl;

  // Momentum sum rule
  std::cout << "\nMomentum sum rule: " << Dists.at(0).Integrate(xmin, 1) + Dists.at(1).Integrate(xmin, 1) << std::endl;

  // Total valence sum rule
  std::cout << "Total valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * Dists.at(2)).Integrate(xmin, 1) << std::endl;

  // V3 valence sum rule
  std::cout << "V3 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * Dists.at(4)).Integrate(xmin, 1) << std::endl;

  // V8 valence sum rule
  std::cout << "V8 valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * Dists.at(6)).Integrate(xmin, 1) << std::endl;

  // Rotate PDFs into the physical basis
  const std::map<int, apfel::Distribution> DistsPhys = apfel::QCDEvToPhys(Dists.GetObjects());

  // Up valence sum rule
  std::cout << "Up valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(2) - DistsPhys.at(-2) )).Integrate(xmin, 1) << std::endl;

  // Down valence sum rule
  std::cout << "Down valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(1) - DistsPhys.at(-1) )).Integrate(xmin, 1) << std::endl;

  // Strange valence sum rule
  std::cout << "Strange valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(3) - DistsPhys.at(-3) )).Integrate(xmin, 1) << std::endl;

  // Get splitting functions and matching functions for the first
  // number of active flavours available
  const std::map<int, apfel::Set<apfel::Operator>> SplitFuncs = DglapObj.begin()->second.SplittingFunctions;
  const std::map<int, apfel::Set<apfel::Operator>> MatchFuncs = DglapObj.begin()->second.MatchingConditions;
  const int nf = DglapObj.begin()->first;

  // Functions to convolute with splitting and matching functions to
  // obtain the appropriate integrals
  const apfel::Distribution MSRDist{g, [] (double const& x) -> double { return 1 / x; }};
  const apfel::Distribution VSRDist{g, [] (double const&) -> double { return 1; }};

  std::cout << "\nChecking sum rules at the level of splitting functions..." << std::endl;

  // Run over orders
  for (auto const& sf : SplitFuncs)
    {
      const std::map<int, apfel::Operator> P = sf.second.GetObjects();
      std::cout << "\norder: " << sf.first << std::endl;
      std::cout << "- Momentum sum rules:" << std::endl;
      std::cout << "  * MSR1 = " << xmin * (( P.at(apfel::EvolutionBasisQCD::PQQ) + P.at(apfel::EvolutionBasisQCD::PGQ) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR2 = " << xmin * (( P.at(apfel::EvolutionBasisQCD::PGG) + P.at(apfel::EvolutionBasisQCD::PQG) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "- Valence sum rules:" << std::endl;
      std::cout << "  * VSR1 = " << (P.at(apfel::EvolutionBasisQCD::PNSV) * VSRDist).Evaluate(xmin) << std::endl;
    }

  std::cout << "\nChecking sum rules at the level of matching functions..." << std::endl;

  // Run over orders
  for (auto const& mf : MatchFuncs)
    {
      const std::map<int, apfel::Operator> K = mf.second.GetObjects();
      std::cout << "\norder: " << mf.first << std::endl;
      std::cout << "- Momentum sum rule:" << std::endl;
      std::cout << "  * MSR1 = " << xmin * (( K.at(apfel::MatchingBasisQCD::M1) + K.at(apfel::MatchingBasisQCD::M4) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR2 = " << xmin / ( nf + 1 ) * (( K.at(apfel::MatchingBasisQCD::M2) - K.at(apfel::MatchingBasisQCD::M3)
                                                           + K.at(apfel::MatchingBasisQCD::M5) - K.at(apfel::MatchingBasisQCD::M6) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR3 (intrinsic heavy flavour) = " << xmin / ( nf + 1 ) * (( K.at(apfel::MatchingBasisQCD::M2) + K.at(apfel::MatchingBasisQCD::M5)
                                                                                     + nf * ( K.at(apfel::MatchingBasisQCD::M3) + K.at(apfel::MatchingBasisQCD::M6) ) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "- Valence sum rule: " << std::endl;
      std::cout << "  * VSR1 = " << (K.at(apfel::MatchingBasisQCD::M7) * VSRDist).Evaluate(xmin) << std::endl;
    }

  return 0;
}
