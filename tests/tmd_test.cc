//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

// Main program
int main()
{
  // Vectors of masses and thresholds
  const std::vector<double> Thresholds = {0, 0, 0, sqrt(2), 4.5, 175};
  //const std::vector<double> Thresholds = {0, 0, 0, 0, 0};

  // Running coupling
  apfel::AlphaQCD a{0.35, sqrt(2), Thresholds, apfel::FixedOrderAccuracy::NNLO};
  const apfel::TabulateObject<double> tabas{a, 500, 0.9, 1001, 3};
  const auto Alphas = [=] (double const& mu) -> double{ return tabas.Evaluate(mu); };

  // x-space grid
  const apfel::Grid g{{apfel::SubGrid{100, 1e-5, 3}, apfel::SubGrid{100, 1e-1, 3}, apfel::SubGrid{100, 6e-1, 3}, apfel::SubGrid{80, 8.5e-1, 5}}};

  // Construct the DGLAP objects
  const auto EvolvedPDFs = BuildDglap(InitializeDglapObjectsQCD(g, Thresholds), apfel::LHToyPDFs, sqrt(2), apfel::FixedOrderAccuracy::NNLO, Alphas);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabPDFs{*EvolvedPDFs, 200, 1, 20000, 3};
  const auto CollPDFs = [=] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabPDFs.Evaluate(mu); };

  // Get timer
  apfel::Timer t;

  // Initialize TMD objects
  const auto TmdObj = apfel::InitializeTmdObjects(g, Thresholds);

  // Kinematics
  const double Vs = 13000;
  const double y  = 0;
  const double Q  = apfel::ZMass;

  // Compute 'x1' and 'x2'
  const double x1 = Q * exp(y) / Vs;
  const double x2 = Q / exp(y) / Vs;

  // Build evolved TMD PDFs
  const auto EvTMDPDFs = BuildTmdPDFs(TmdObj, CollPDFs, Alphas, apfel::LogAccuracy::NNNLL);

  // Get Drell-Yan hard-factor function
  const double hcs = HardFactor("DY", TmdObj, Alphas, apfel::LogAccuracy::NNNLL)(Q);

  // Number of active flavours at 'Q'
  const int nf = apfel::NF(Q, Thresholds);

  // EW charges
  const std::vector<double> Bq = apfel::ElectroWeakCharges(Q, true);

  // Get Evolved TMD PDFs at bmax. This is is subtracted off from the
  // b-dependent luminosity to improve the convergence of the Hankel
  // transform, especially at small qT.
  const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(2 * exp( - apfel::emc), Q, Q * Q).GetObjects());
  double lumimax = 0;
  for (int i = 1; i <= nf; i++)
    lumimax += Bq[i-1] * ( xF.at(i).Evaluate(x1) * xF.at(-i).Evaluate(x2) + xF.at(-i).Evaluate(x1) * xF.at(i).Evaluate(x2) );

  // Construct the TMD luminosity in b space to be trasformed in qT
  // space.
  const std::function<double(double const&)> TMDLumib = [=] (double const& b) -> double
  {
    // Compute b*
    const double bstar = std::max(b / sqrt( 1 + pow(b / 2 / exp( - apfel::emc), 2) ), 1e-4);

    // Get Evolved TMD PDFs and rotate them into the physical
    // basis
    const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(bstar, Q, Q * Q).GetObjects());

    // Combine TMDs through the EW charges
    double lumi = 0;
    for (int i = 1; i <= nf; i++)
      lumi += Bq[i-1] * ( xF.at(i).Evaluate(x1) * xF.at(-i).Evaluate(x2) + xF.at(-i).Evaluate(x1) * xF.at(i).Evaluate(x2) );

    // Combine all pieces and return
    return b * ( lumi - lumimax );
  };

  // Double exponential quadrature
  //apfel::DoubleExponentialQuadrature DEObj{};
  apfel::OgataQuadrature DEObj{};

  // Phase-space reduction factor
  apfel::TwoBodyPhaseSpace ps{20, -1, 2.4};

  // Compute predictions
  const int nqT = 100;
  const double qTmin = 0.01;
  const double qTmax = 30;
  const double qTstp = ( qTmax - qTmin ) / ( nqT - 1 );
  for (double qT = qTmin; qT <= qTmax; qT += qTstp)
    std::cout << std::scientific << Q << "  " << y << "  " << qT << "  "
              << apfel::ConvFact * qT * 8 * M_PI * pow(apfel::alphaem, 2) * hcs / pow(Q, 3) / 9 * DEObj.transform(TMDLumib, qT) << "  "
              << ps.PhaseSpaceReduction(Q, y, qT) << "  "
              << ps.ParityViolatingPhaseSpaceReduction(Q, y, qT) << "  "
              << std::endl;
  t.stop();

  return 0;
}
