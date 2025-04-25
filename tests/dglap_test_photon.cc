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
  const auto as  = [&] (double const& mu) -> double{ return Alphas.Evaluate(mu); };
  const auto aem = [&] (double const& mu) -> double{ return Alpha.Evaluate(mu); };

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
  const std::map<int, apfel::Distribution> oppdfs = apfel::PlusMinusQCDQEDToPhys((tops * pdfs0).GetObjects());

  // Get PDFs at the final scale as distributions
  const apfel::Set<apfel::Distribution> Dists = TabulatedPDFs.Evaluate(mu);

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

  std::cout << "Direct Evolution:" << std::endl;
  for (auto const& x : xlha)
    {
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << (pdfs.at(2) - pdfs.at(-2)).Evaluate(x)
                << "  " << (pdfs.at(1) - pdfs.at(-1)).Evaluate(x)
                << "  " << 2 * (pdfs.at(-2) + pdfs.at(-1)).Evaluate(x)
                << "  " << (pdfs.at(4) + pdfs.at(-4)).Evaluate(x)
                << "  " << pdfs.at(0).Evaluate(x)
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

  std::cout << "Evolution through the interpolated evolution operator:" << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> opxpdfs = apfel::PlusMinusQCDQEDToPhys((apfel::Set<apfel::Distribution> {apfel::EvolveDistributionsBasisQCDQED{}, tops.Evaluate(x).GetObjects()} * pdfs0).Squash());
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << opxpdfs.at(2) - opxpdfs.at(-2)
                << "  " << opxpdfs.at(1) - opxpdfs.at(-1)
                << "  " << 2 * ( opxpdfs.at(-2) + opxpdfs.at(-1) )
                << "  " << opxpdfs.at(4) + opxpdfs.at(-4)
                << "  " << opxpdfs.at(0)
                << std::endl;
    }
  std::cout << "\n";

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

  std::cout << "Interpolation on the PDF table as a map (x and Q independently):" << std::endl;
  for (auto const& x : xlha)
    {
      const std::map<int, double> DistMap = apfel::PlusMinusQCDQEDToPhys(TabulatedPDFs.EvaluateMapxQ(x, mu));
      std::cout.precision(1);
      std::cout << x;
      std::cout.precision(4);
      std::cout << "  " << DistMap.at(2) - DistMap.at(-2)
                << "  " << DistMap.at(1) - DistMap.at(-1)
                << "  " << 2 * ( DistMap.at(-2) + DistMap.at(-1) )
                << "  " << DistMap.at(4) + DistMap.at(-4)
                << "  " << DistMap.at(0)
                << std::endl;
    }
  std::cout << "\n";

  int k = 1000000;
  std::cout << "Interpolating " << k << " times a single PDF on the (x,Q) grid... ";
  t.start();
  for (int i = 0; i < k; i++)
    TabulatedPDFs.EvaluatexQ(0, 0.05, mu);
  t.stop();

  k = 100000;
  std::cout << "Interpolating " << k << " times a map of PDFs on the (x,Q) grid... ";
  t.start();
  for (int i = 0; i < k; i++)
    TabulatedPDFs.EvaluateMapxQ(0.05, mu);
  t.stop();

  return 0;
}
