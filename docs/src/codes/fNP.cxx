//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/apfelxx.h>
#include <NangaParbat/bstar.h>

#include "LHAPDF/LHAPDF.h"

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
      distpdf->quarkThreshold(5),
      1e3};

  // Hard scale
  const double Q = 91;

  // Alpha_s.
  //const auto Alphas = [&] (double const& mu) -> double{ return distpdf->alphasQ(mu); };
  apfel::AlphaQCD a{distpdf->alphasQ(Q), Q, Thresholds, 2};
  const apfel::TabulateObject<double> TAlphas{a, 100, 0.49999, 1001, 3};
  const auto Alphas = [&] (double const& mu) -> double{ return TAlphas.Evaluate(mu); };

  // Perturbative order.
  const int PerturbativeOrder = apfel::LogAccuracy::NNNLL;

  // Initial scale-variation factor;
  const double Ci = 1;

  // Final scale-variation factor;
  const double Cf = 1;

  // x-space grid.
  const apfel::Grid g{{{100,1e-5,3}, {60,1e-1,3}, {50,6e-1,3}, {100,8e-1,3}}};

  // Initialize QCD evolution objects
  const auto DglapObj = InitializeDglapObjectsQCD(g, Thresholds);

  // Rotate PDF set into the QCD evolution basis.
  const auto RotPDFs = [=] (double const& x, double const& mu) -> std::map<int,double>{ return apfel::PhysToQCDEv(distpdf->xfxQ(x,mu)); };

  // Construct the DGLAP objects
  const auto EvolvedPDFs2 = BuildDglap(DglapObj, RotPDFs, Q, 2, Alphas);

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabulatedPDFs{*EvolvedPDFs2, 100, 0.5, 1000, 3};

  // Construct set of distributions as a function of the scale to be
  // tabulated.
  const auto EvolvedPDFs = [=,&g] (double const& mu) -> apfel::Set<apfel::Distribution>{
    return apfel::Set<apfel::Distribution>{apfel::EvolutionBasisQCD{apfel::NF(mu, Thresholds)}, DistributionMap(g, RotPDFs, mu)};
  };

  // Tabulate PDFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabPDFs{EvolvedPDFs, 50, 1, 100000, 3, Thresholds};
  //const auto CollPDFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabPDFs.Evaluate(mu); };
  const auto CollPDFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabulatedPDFs.Evaluate(mu); };

  // Initialize TMD objects.
  const auto TmdObj = InitializeTmdObjects(g, Thresholds);

  // Get evolution factors.
  const auto QuarkEvolFactor = QuarkEvolutionFactor(TmdObj, Alphas, PerturbativeOrder, Ci);

  // Match TMD PDFs
  const auto MatchedTMDPDFs = MatchTmdPDFs(TmdObj, CollPDFs, Alphas, PerturbativeOrder, Ci);

  // TMD scales
  const double muf   = Cf * Q;
  const double zetaf = Q * Q;

  // Construct function that returns evolved TMDs including the
  // non-perturbative part. This can be tabulated in b.
  const auto TMD = [=] (double const& b) -> apfel::Set<apfel::Distribution>{ return QuarkEvolFactor(b, muf, zetaf) * MatchedTMDPDFs(b); };

  const double xb = 0.98;

  const int nb = 50;
  const double bmin = 0.1;
  const double bmax = 1;
  const double bstep = ( bmax - bmin ) / ( nb - 1 );

  const double bTmin  = 2 * exp( - apfel::emc) / sqrt(zetaf);
  const double power = 4;

  double b = bmin;
  std::cout << std::scientific;
  for (int ib = 0; ib < nb; ib++)
    {
      const double bc = b / pow( ( 1 - exp( - pow(b / bTmin, power) ) ), 1 / power);
      std::cout << 2 * exp( - apfel::emc) / b << "\t"
		<< b << "\t" << TMD(bc).at(1).Evaluate(xb) / TMD(NangaParbat::bstarmin(b, Q)).at(1).Evaluate(xb) << "\t"
		<< TMD(b).at(1).Evaluate(xb) << "\t"
		<< TMD(NangaParbat::bstarmin(b, Q)).at(1).Evaluate(xb) << std::endl;
      b += bstep;
    }
  return 0;
}
