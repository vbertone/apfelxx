//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

int main()
{
  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Initial scale
  const double mu0 = sqrt(2);

  // Vector of thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};

  // Perturbative order
  const int PerturbativeOrder = 2;

  // Running strong and electromagnetic couplings
  apfel::AlphaQCD a0{0.35, sqrt(2), Thresholds, PerturbativeOrder};
  apfel::AlphaQED a1{7.496252e-3, sqrt(2), Thresholds, {}, 1};
  const apfel::TabulateObject<double> Alphas{a0, 100, 0.9, 1001, 3};
  const apfel::TabulateObject<double> Alpha{a1, 100, 0.9, 1001, 3};
  const auto as  = [=] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
  const auto aem = [=] (double const& mu) -> double{ return Alpha.Evaluate(mu); };

  // Initialize QCD evolution objects
  const auto DglapObj   = InitializeDglapObjectsPhoton(g, Thresholds);
  const auto DglapObjOp = InitializeDglapObjectsPhoton(g, Thresholds, true);

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(DglapObj, apfel::LHToyPDFsQCDQED, mu0, PerturbativeOrder, as, aem);
  const auto EvolvedOps  = BuildDglap(DglapObjOp,                       mu0, PerturbativeOrder, as, aem);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs, 50, 1, 1000, 3};

  // Tabulate Operators
  const apfel::TabulateObject<apfel::Set<apfel::Operator>> TabulatedOps{*EvolvedOps, 50, 1, 1000, 3};

  // Final scale
  const double mu = 100;

  // Compute results
  std::cout << "Direct evolution (4th order Runge-Kutta) from Q0 = " << mu0 << " GeV to Q = " << mu << " GeV... ";

  // Evolve PDFs to the final Scale
  apfel::Timer t;
  const std::map<int, apfel::Distribution> pdfs = apfel::PlusMinusQCDQEDToPhys(EvolvedPDFs->Evaluate(mu).GetObjects());
  t.stop();

  std::cout << "Interpolation of the tabulated PDFs... ";
  t.start();
  const std::map<int, apfel::Distribution> tpdfs = apfel::PlusMinusQCDQEDToPhys(TabulatedPDFs.Evaluate(mu).GetObjects());
  t.stop();

  std::cout << "Interpolation of the tabulated evolution operators... ";
  t.start();
  apfel::Set<apfel::Operator> tops = TabulatedOps.Evaluate(mu);
  t.stop();

  // Set appropriate convolution basis for the evolution operators and
  // convolute them with initial-scale distributions.
  tops.SetMap(apfel::EvolveDistributionsBasisQCDQED{});
  const apfel::Set<apfel::Distribution> pdfs0{apfel::EvolveDistributionsBasisQCDQED{}, DistributionMap(g, apfel::LHToyPDFsQCDQED, mu0)};

  // Compute inhomogeneous term to include when evolving PDFs through
  // the evolution operator.
  const auto InHom = InhomogeneousTermsQCDQED(DglapObj, PerturbativeOrder, as, aem);
  const std::function<apfel::Set<apfel::Distribution>(double const&)> GammaInHom = [=] (double const& mup) -> apfel::Set<apfel::Distribution>
  {
    const int nf = apfel::NF(mup, Thresholds);
    apfel::Set<apfel::Operator> gamma = BuildDglap(DglapObjOp, mup, PerturbativeOrder, as, aem)->Evaluate(mu);
    apfel::Set<apfel::Distribution> inhom = InHom(nf, 2 * log(mup));
    gamma.SetMap(apfel::EvolveDistributionsBasisQCDQED{});
    inhom.SetMap(apfel::EvolveDistributionsBasisQCDQED{});
    return ( 2 / mup ) * gamma * inhom;
  };
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabGammaInHom{GammaInHom, 100, mu0, mu, 3, Thresholds};

  apfel::Set<apfel::Distribution> IntGammaInHom = TabGammaInHom.Integrate(mu0, mu);

  IntGammaInHom.SetMap(apfel::EvolveDistributionsBasisQCDQED{});
  const std::map<int, apfel::Distribution> oppdfs = apfel::PlusMinusQCDQEDToPhys((tops * pdfs0 + IntGammaInHom).GetObjects());

  // Print results
  const std::vector<double> xlha = {1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 3e-1, 5e-1, 7e-1, 9e-1};
  std::cout << std::scientific;

  std::cout << "\nAlphaQCD(Q) = " << as(mu) << std::endl;
  std::cout << "AlphaQED(Q) = " << aem(mu) << std::endl;
  std::cout << "\n   x    "
            << "   u-ubar   "
            << "   d-dbar   "
            << " 2(ubr+dbr) "
            << "   c+cbar   "
            << "    gluon   "
            << std::endl;

  std::cout << "Interpolation on the PDF table (all x for each Q):" << std::endl;
  for (auto const& x : xlha)
    {
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << (tpdfs.at(2) - tpdfs.at(-2)).Evaluate(x)
                << "  " << (tpdfs.at(1) - tpdfs.at(-1)).Evaluate(x)
                << "  " << 2 * (tpdfs.at(-2) + tpdfs.at(-1)).Evaluate(x)
                << "  " << (tpdfs.at(4) + tpdfs.at(-4)).Evaluate(x)
                << "  " << tpdfs.at(0).Evaluate(x)
                << std::endl;
    }
  std::cout << "\n";

  std::cout << "Evolution through the evolution operator:" << std::endl;
  for (auto const& x : xlha)
    {
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << (oppdfs.at(2) - oppdfs.at(-2)).Evaluate(x)
                << "  " << (oppdfs.at(1) - oppdfs.at(-1)).Evaluate(x)
                << "  " << 2 * (oppdfs.at(-2) + oppdfs.at(-1)).Evaluate(x)
                << "  " << (oppdfs.at(4) + oppdfs.at(-4)).Evaluate(x)
                << "  " << oppdfs.at(0).Evaluate(x)
                << std::endl;
    }
  std::cout << "\n";

  return 0;
}
