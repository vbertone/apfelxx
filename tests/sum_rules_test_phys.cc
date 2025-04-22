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
  const auto DglapObj = InitializeDglapObjectsQCDPhys(g, Thresholds, false, apfel::eps5, true, {0, 0, 0, 0, 0, 0, 0});

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFsPlusMinus, mu0, PerturbativeOrder, as);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Get PDFs at the final scale as distributions
  const apfel::Set<apfel::Distribution> Dists = TabulatedPDFs.Evaluate(mu);

  // Lower integration limit
  const double xmin = 1e-9;

  // Start the checks
  std::cout << "\nChecking sum rules at the level of PDFs..." << std::endl;

  // Momentum sum rule
  std::cout << "\nMomentum sum rule: " << (Dists.at(6) + Dists.at(7) + Dists.at(8) + Dists.at(9) + Dists.at(10) + Dists.at(11) + Dists.at(12)).Integrate(xmin, 1) << std::endl;

  // Total valence sum rule
  std::cout << "Total valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( Dists.at(0) + Dists.at(1) + Dists.at(2) + Dists.at(3) + Dists.at(4) + Dists.at(5) )).Integrate(xmin, 1) << std::endl;

  // Rotate PDFs into the physical basis
  const std::map<int, apfel::Distribution> DistsPhys = apfel::PlusMinusToPhys(Dists.GetObjects());

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
      std::cout << "  * MSR1 = " << xmin * (( P.at(apfel::PhysicalBasisQCD::PNS) + ( nf - 1 ) * P.at(apfel::PhysicalBasisQCD::PPS) + P.at(apfel::PhysicalBasisQCD::PGQ) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR2 = " << xmin * (( P.at(apfel::PhysicalBasisQCD::PGG) + nf * P.at(apfel::PhysicalBasisQCD::PQG) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "- Valence sum rules:" << std::endl;
      std::cout << "  * VSR1 = " << (( P.at(apfel::PhysicalBasisQCD::PNV) + ( nf - 1 ) * P.at(apfel::PhysicalBasisQCD::PPV) ) * VSRDist).Evaluate(xmin) << std::endl;
    }

  std::cout << "\nChecking sum rules at the level of matching functions..." << std::endl;

  // Run over orders
  for (auto const& mf : MatchFuncs)
    {
      if (mf.first < 0)
        continue;
      const std::map<int, apfel::Operator> K = mf.second.GetObjects();
      std::cout << "\norder: " << mf.first << std::endl;
      std::cout << "- Momentum sum rule:" << std::endl;
      std::cout << "  * MSR1 = " << xmin * (( K.at(apfel::PhysicalMatchingBasisQCD::KGG) + nf * K.at(apfel::PhysicalMatchingBasisQCD::KLG) + K.at(apfel::PhysicalMatchingBasisQCD::KHG) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR2 = " << xmin * (( K.at(apfel::PhysicalMatchingBasisQCD::KGL) + K.at(apfel::PhysicalMatchingBasisQCD::KLL)
                                              + nf * K.at(apfel::PhysicalMatchingBasisQCD::KLLP) + K.at(apfel::PhysicalMatchingBasisQCD::KHL) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR3 (intrinsic heavy flavour) = " << xmin * (( K.at(apfel::PhysicalMatchingBasisQCD::KGH) + K.at(apfel::PhysicalMatchingBasisQCD::KHH) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "- Valence sum rule: " << std::endl;
      std::cout << "  * VSR1 = " << (K.at(apfel::PhysicalMatchingBasisQCD::KLL) * VSRDist).Evaluate(xmin) << std::endl;
    }

  return 0;
}
