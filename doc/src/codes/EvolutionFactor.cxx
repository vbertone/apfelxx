//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>

#include "LHAPDF/LHAPDF.h"

// Non-perturnative function
double fNP(double const& b, double const& zetaf)
{
  const double g1 = 0.02;
  const double g2 = 0.5;
  const double Q0 = 1.6;
  return exp( ( - g1 - g2 * log( sqrt(zetaf) / Q0 / 2 ) ) * b * b / 2 );
}

// b* prescription
double bstar(double const& b)
{
  const double bmax = 1;
  const double bs = b / sqrt( 1 + pow(b / bmax, 2) );
  return bs;
}

// Main program
int main()
{
 // Open LHAPDF set.
  LHAPDF::PDF* distpdf = LHAPDF::mkPDF("MMHT2014nnlo68cl");

  // Thresholds.
  const std::vector<double> Thresholds{distpdf->quarkThreshold(1),
      distpdf->quarkThreshold(2),
      distpdf->quarkThreshold(3),
      distpdf->quarkThreshold(4),
      distpdf->quarkThreshold(5)};

  // Alpha_s.
  const auto Alphas = [&] (double const& mu) -> double{ return distpdf->alphasQ(mu); };

  // Perturbative order.
  const int PerturbativeOrder = apfel::LogAccuracy::NNLL;

  // Initial scale-variation factor;
  const double Ci = 1;

  // Final scale-variation factor;
  const double Cf = 1;

  // x-space grid.
  const apfel::Grid g{{{100,1e-5,3}, {60,1e-1,3}, {50,6e-1,3}, {100,8e-1,3}}};

  // Rotate PDF set into the QCD evolution basis.
  const auto RotPDFs = [=] (double const& x, double const& mu) -> std::map<int,double>{ return apfel::PhysToQCDEv(distpdf->xfxQ(x,mu)); };

  // Construct set of distributions as a function of the scale to be
  // tabulated.
  const auto EvolvedPDFs = [=,&g] (double const& mu) -> apfel::Set<apfel::Distribution>{
    return apfel::Set<apfel::Distribution>{apfel::EvolutionBasisQCD{apfel::NF(mu, Thresholds)}, DistributionMap(g, RotPDFs, mu)};
  };

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabPDFs{EvolvedPDFs, 50, 1, 100000, 3, Thresholds};
  const auto CollPDFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabPDFs.Evaluate(mu); };

  // Initialize TMD objects.
  const auto TmdObj = InitializeTmdObjects(g, Thresholds);

  // Get evolution factors.
  const auto QuarkEvolFactor = QuarkEvolutionFactor(TmdObj, Alphas, PerturbativeOrder, Ci);

  // Match TMD PDFs
  const auto MatchedTMDPDFs = MatchTmdPDFs(TmdObj, CollPDFs, Alphas, PerturbativeOrder, Ci);

  // Hard scale
  const double Q = 91.1876;

  // Rapidity
  const double y = 0;

  // C.M.E.
  const double Vs = 13000;

  // TMD scales
  const double muf   = Cf * Q;
  const double zetaf = Q * Q;

  // Construct function that returns evolved TMDs including the
  // non-perturbative part. This can be tabulated in b.
  const auto EvolvedTMDPDFs = [=] (double const& b) -> apfel::Set<apfel::Distribution>
    {
      const double bs = bstar(b);
      return fNP(b, zetaf) * QuarkEvolFactor(bs, muf, zetaf) * MatchedTMDPDFs(bs);
    };

  // Tabulate input TMDs in the impact parameter to make the
  // integral faster.
  const auto TabFunc    = [] (double const& b) -> double{ return log(b); };
  const auto InvTabFunc = [] (double const& fb) -> double{ return exp(fb); };
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabEvolvedTMDPDFs{EvolvedTMDPDFs, 50, 1e-4, 10, 3, {}, TabFunc, InvTabFunc};

  // EW charges
  const std::vector<double> Bq = apfel::ElectroWeakCharges(Q, true);

  // Compute hard coefficient, including the other kinematic
  // factors.
  const double alpha = 1. / 127.;
  const double hcs = apfel::ConvFact * 8 * M_PI * pow(alpha, 2) * apfel::HardFactorDY(PerturbativeOrder, Alphas(muf), apfel::NF(muf, Thresholds), Cf) / 9 / pow(Q, 3);

  // Ogata quadrature object
  apfel::OgataQuadrature bintegrand{};

  // Define qT-distribution function 
  const auto qTdist = [=] (double const& qT) -> double
  {
    // Bjorken variables
    const double x1 = Q * exp(-y) / Vs;
    const double x2 = Q * exp(y) / Vs;

    // Construct the TMD luminosity in b scale to be fed to be
    // trasformed in qT space.
    const auto TMDLumib = [=] (double const& b) -> double
    {
      // Get map of the TMDs in "x1" and "x2" and rotate them into the
      // physical basis.
      std::map<int,double> TMD1 = apfel::QCDEvToPhys(TabEvolvedTMDPDFs.EvaluateMapxQ(x1, b));
      std::map<int,double> TMD2 = apfel::QCDEvToPhys(TabEvolvedTMDPDFs.EvaluateMapxQ(x2, b));

      // Construct the combination of TMDs weighted by the EW
      // charges. Remember that each TMD has a factor x in the
      // front. This effectively means that the luminosity is
      // multiplied by x1 * x2 = Q2 / s.
      double lumi = 0;
      for (int i = 1; i <= 5; i++)
	lumi += Bq[i-1] * ( TMD1[i] * TMD2[i] + TMD1[-i] * TMD2[-i] );

	// Multiply by "b" and divide by two to reduce to the Fourier
	// transform to a Hankel transform.
	return b * lumi / 2;
    };
    return 2 * qT * hcs * bintegrand.transform(TMDLumib, qT);
  };

  apfel::Timer t;
  const std::vector<double> qTv{0.2500e+00, 0.7500e+00, 0.1250e+01, 0.1750e+01, 0.2250e+01, 0.2750e+01, 0.3250e+01, 0.3750e+01,
      0.4250e+01, 0.4750e+01, 0.5250e+01, 0.5750e+01, 0.6250e+01, 0.6750e+01, 0.7250e+01, 0.7750e+01, 0.8250e+01,
      0.8750e+01, 0.9250e+01, 0.9750e+01, 0.1025e+02, 0.1075e+02, 0.1125e+02, 0.1175e+02, 0.1250e+02, 0.1350e+02,
      0.1450e+02, 0.1550e+02, 0.1650e+02, 0.1750e+02, 0.1850e+02};
  std::cout << std::scientific << "  qT [GeV]    "
       << "   sigma      "
       << std::endl;
  for (auto const& qT : qTv)
    std::cout << qT << "  " << qTdist(qT) << std::endl;
  std::cout << "\n";
  t.stop();

  return 0;
}
