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

  // Vector of thresholds
  const std::vector<double> QuarkThresholds = {0, 0, 0, sqrt(2), 4.5, 175};
  const std::vector<double> LeptonThresholds = {0, 0, 1.777};

  // Perturbative order
  const int PerturbativeOrder = 3;

  // Running couplings
  apfel::AlphaQCDQED a{0.35, 7.496252e-3, sqrt(2), QuarkThresholds, LeptonThresholds, PerturbativeOrder};
  const apfel::TabulateObject<apfel::matrix<double>> Couplings{a, 100, 0.9, 1001, 3};
  const auto as  = [&] (double const& mu) -> double{ return Couplings.Evaluate(mu)(0, 0); };
  const auto aem = [&] (double const& mu) -> double{ return Couplings.Evaluate(mu)(1, 0); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{200, 1e-9, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Initialise DGLAP objects
  const auto DglapObj = InitializeDglapObjectsQCDQED(g, QuarkThresholds, LeptonThresholds, false, apfel::eps5, true, {0, 0, 0, 0});

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFsQCDQED, mu0, PerturbativeOrder, as, aem);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Get PDFs at the final scale as distributions
  const apfel::Set<apfel::Distribution> Dists = TabulatedPDFs.Evaluate(mu);

  // Lower integration limit
  const double xmin = 1e-9;

  // Start the checks
  std::cout << "\nChecking sum rules at the level of PDFs..." << std::endl;

  // Momentum sum rule
  std::cout << "\nMomentum sum rule: " << (Dists.at(0) + Dists.at(1) + Dists.at(2) + Dists.at(3) + Dists.at(4) + Dists.at(5) + Dists.at(6) + Dists.at(7) + Dists.at(8) + Dists.at(9) + Dists.at(10)).Integrate(xmin, 1) << std::endl;

  // Total valence sum rule
  std::cout << "Total valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( Dists.at(11) + Dists.at(12) + Dists.at(13) + Dists.at(14) + Dists.at(15) + Dists.at(16) + Dists.at(17) + Dists.at(18) + Dists.at(19) )).Integrate(xmin, 1) << std::endl;

  // Rotate PDFs into the physical basis
  const std::map<int, apfel::Distribution> DistsPhys = apfel::PlusMinusQCDQEDToPhys(Dists.GetObjects());

  // Up valence sum rule
  std::cout << "Up valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(2) - DistsPhys.at(-2) )).Integrate(xmin, 1) << std::endl;

  // Down valence sum rule
  std::cout << "Down valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(1) - DistsPhys.at(-1) )).Integrate(xmin, 1) << std::endl;

  // Strange valence sum rule
  std::cout << "Strange valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(3) - DistsPhys.at(-3) )).Integrate(xmin, 1) << std::endl;

  // Electron valence sum rule
  std::cout << "Electron valence sum rule: " << ([] (double const& x) -> double { return 1 / x; } * ( DistsPhys.at(11) - DistsPhys.at(-11) )).Integrate(xmin, 1) << std::endl;

  // Get splitting functions and matching functions for the first
  // number of active flavours available
  const std::map<std::pair<int, int>, apfel::Set<apfel::Operator>> SplitFuncs = DglapObj.begin()->second.SplittingFunctions;
  const std::map<std::pair<int, int>, apfel::Set<apfel::Operator>> MatchFuncs = DglapObj.begin()->second.MatchingConditions;
  const std::vector<int> actfl = DglapObj.begin()->second.ActiveFlavours;
  //const int nt = DglapObj.begin()->first;

  // Functions to convolute with splitting and matching functions to
  // obtain the appropriate integrals
  const apfel::Distribution MSRDist{g, [] (double const& x) -> double { return 1 / x; }};
  const apfel::Distribution VSRDist{g, [] (double const&) -> double { return 1; }};

  std::cout << "\nChecking sum rules at the level of splitting functions..." << std::endl;

  // Run over orders
  for (auto const& sf : SplitFuncs)
    {
      std::cout << "\norder: [" << sf.first.first << ", " << sf.first.second << "]" << std::endl;
      const int nD = actfl[0];
      const int nU = actfl[1];
      const int nL = actfl[2];
      const std::map<int, apfel::Operator> P = sf.second.GetObjects();
      std::cout << "- Momentum sum rules:" << std::endl;
      std::cout << "  * MSR1 = " << xmin * (( P.at(apfel::EvolutionBasisQCDQED::PPLL) + nL * P.at(apfel::EvolutionBasisQCDQED::PPSLL) + nU * P.at(apfel::EvolutionBasisQCDQED::PPSUL) + nD * P.at(apfel::EvolutionBasisQCDQED::PPSDL) + P.at(apfel::EvolutionBasisQCDQED::PgmL) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR2 = " << xmin * (( P.at(apfel::EvolutionBasisQCDQED::PPUU) + nL * P.at(apfel::EvolutionBasisQCDQED::PPSLU) + nU * P.at(apfel::EvolutionBasisQCDQED::PPSUU) + nD * P.at(apfel::EvolutionBasisQCDQED::PPSDU) + P.at(apfel::EvolutionBasisQCDQED::PgU) + P.at(apfel::EvolutionBasisQCDQED::PgmU) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR3 = " << xmin * (( P.at(apfel::EvolutionBasisQCDQED::PPDD) + nL * P.at(apfel::EvolutionBasisQCDQED::PPSLD) + nU * P.at(apfel::EvolutionBasisQCDQED::PPSUD) + nD * P.at(apfel::EvolutionBasisQCDQED::PPSDD) + P.at(apfel::EvolutionBasisQCDQED::PgD) + P.at(apfel::EvolutionBasisQCDQED::PgmD) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR4 = " << xmin * (( nU * P.at(apfel::EvolutionBasisQCDQED::PUg) + nD * P.at(apfel::EvolutionBasisQCDQED::PDg) + P.at(apfel::EvolutionBasisQCDQED::Pgmg) + P.at(apfel::EvolutionBasisQCDQED::Pgg) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR5 = " << xmin * (( nL * P.at(apfel::EvolutionBasisQCDQED::PLgm) + nU * P.at(apfel::EvolutionBasisQCDQED::PUgm) + nD * P.at(apfel::EvolutionBasisQCDQED::PDgm) + P.at(apfel::EvolutionBasisQCDQED::Pggm) + P.at(apfel::EvolutionBasisQCDQED::Pgmgm) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "- Valence sum rules:" << std::endl;
      std::cout << "  * VSR1 = " << (( P.at(apfel::EvolutionBasisQCDQED::PMDD) + nD * P.at(apfel::EvolutionBasisQCDQED::PPV) ) * VSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * VSR2 = " << (( P.at(apfel::EvolutionBasisQCDQED::PMUU) + nU * P.at(apfel::EvolutionBasisQCDQED::PPV) ) * VSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * VSR3 = " << (P.at(apfel::EvolutionBasisQCDQED::PMLL) * VSRDist).Evaluate(xmin) << std::endl;
    }

  std::cout << "\nChecking sum rules at the level of matching functions..." << std::endl;

  // Run over orders
  for (auto const& mf : MatchFuncs)
    {
      if (mf.first.first < 0)
        continue;
      std::cout << "\norder: [" << mf.first.first << ", " << mf.first.second << "]" << std::endl;
      const int nf = actfl[0] + actfl[1];
      const std::map<int, apfel::Operator> K = mf.second.GetObjects();
      std::cout << "- Momentum sum rule:" << std::endl;
      std::cout << "  * MSR1 = " << xmin * (( K.at(apfel::MatchingBasisQCDQED::KQg) + nf * K.at(apfel::MatchingBasisQCDQED::Kqg) + K.at(apfel::MatchingBasisQCDQED::Kgg) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR2 = " << xmin * (( K.at(apfel::MatchingBasisQCDQED::KXgm) + K.at(apfel::MatchingBasisQCDQED::KXgmgm) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR3 = " << xmin * (( K.at(apfel::MatchingBasisQCDQED::KQqp) + K.at(apfel::MatchingBasisQCDQED::KNSq) + nf * K.at(apfel::MatchingBasisQCDQED::Kqqp) + K.at(apfel::MatchingBasisQCDQED::Kgq) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * MSR4 (intrinsic heavy flavour) = " << xmin * (( K.at(apfel::MatchingBasisQCDQED::KXX) + K.at(apfel::MatchingBasisQCDQED::KgQ) + K.at(apfel::MatchingBasisQCDQED::KgmX) ) * MSRDist).Evaluate(xmin) << std::endl;
      std::cout << "- Valence sum rule: " << std::endl;
      std::cout << "  * VSR1 = " << (K.at(apfel::MatchingBasisQCDQED::KNSqm) * VSRDist).Evaluate(xmin) << std::endl;
      std::cout << "  * VSR2 = " << (K.at(apfel::MatchingBasisQCDQED::KNSsqm) * VSRDist).Evaluate(xmin) << std::endl;
    }

  return 0;
}
