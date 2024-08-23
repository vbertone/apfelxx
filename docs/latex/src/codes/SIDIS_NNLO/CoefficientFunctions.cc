#include "CoefficientFunctions.h"
#include <apfel/SIDIS.h>
#include <LHAPDF/LHAPDF.h>

int main()
{
  // PDF and FF sets
  LHAPDF::PDF* pdfs = LHAPDF::mkPDF("MSHT20nnlo_as118");
  LHAPDF::PDF* ffs = LHAPDF::mkPDF("MAPFF10NNLOPIp");

  // Get thresholds from PDF set
  std::vector<double> Thresholds;
  for (auto const& v : pdfs->flavors())
    if (v > 0 && v < 7)
      Thresholds.push_back(pdfs->quarkThreshold(v));

  // Structure functions using APFEL++
  // Define grid
  const apfel::Grid g{{{200, 7e-3, 3}, {180, 2e-1, 3}, {150, 5e-1, 3}, {50, 8.5e-1, 3}}};

  // Initialize SIDIS objects
  const apfel::SidisObjects so = InitializeSIDIS(g, Thresholds);

  // Initialize inclusive structure functions.
  const auto InPDFs = [&] (double const& x, double const& mu) -> std::map<int, double>{ return pdfs->xfxQ(x, mu); };
  const auto InFFs  = [&] (double const& z, double const& mu) -> std::map<int, double>{ return ffs->xfxQ(z, mu); };

  // Define function to compute SIDIS cross section at O(as)
  const std::function<apfel::DoubleObject<apfel::Distribution>(double const&)> CrossSectionFO = [=, &g] (double const& Q) -> apfel::DoubleObject<apfel::Distribution>
    {
      // Compute number of active flavours
      const int nf = apfel::NF(Q, Thresholds);

      // Coupling from PDF set
      const double coup = pdfs->alphasQ(Q) / apfel::FourPi;

      // Compute distributions for PDFs and FFs
      const std::map<int, apfel::Distribution> dPDF = apfel::DistributionMap(g, InPDFs, Q);
      const std::map<int, apfel::Distribution> dFF  = apfel::DistributionMap(g, InFFs,  Q);

      apfel::DoubleObject<apfel::Distribution> distqq;
      apfel::DoubleObject<apfel::Distribution> distgq;
      apfel::DoubleObject<apfel::Distribution> distqg;
      for (auto j = - nf; j <= nf; j++)
	{
	  // Skip the gluon
	  if (j == 0)
	    continue;

	  distqq.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(j),  dFF.at(j)});
	  distgq.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(j),  dFF.at(21)});
	  distqg.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(21), dFF.at(j)});
	}

      // Assemple double distribution for the reduced cross section as Y^+ F2 - y^2 FL
      //return ( so.C20qq + coup * so.C21qq ) * distqq + coup * ( so.C21gq * distgq + so.C21qg * distqg );
      return so.C20qq * distqq;
    };
  const apfel::TabulateObject<apfel::DoubleObject<apfel::Distribution>> TabCrossSectionFO{CrossSectionFO, 50, 1, 10, 3, Thresholds};

  // Construct parton luminosities
  const std::function<double(double const&, double const&, double const&)> LqqNS = [=] (double const& x, double const& z, double const& mu) -> double
    {
      double lqqns = 0;
      for (int i = 1; i <= apfel::NF(mu, Thresholds); i++)
	lqqns += apfel::QCh2[i-1] * ( ffs->xfxQ(i, z, mu) * pdfs->xfxQ(i, x, mu) + ffs->xfxQ(-i, z, mu) * pdfs->xfxQ(-i, x, mu) );
      return lqqns;
    };
  const std::function<double(double const&, double const&, double const&)> LqqPS = [=] (double const& x, double const& z, double const& mu) -> double
    {
      const int Nf = apfel::NF(mu, Thresholds);
      double lqqps = 0;
      for (int i = 1; i <= Nf; i++)
	lqqps += ffs->xfxQ(i, z, mu) * pdfs->xfxQ(i, x, mu) + ffs->xfxQ(-i, z, mu) * pdfs->xfxQ(-i, x, mu);
      return apfel::SumCh2[Nf] * lqqps;
    };
  const std::function<double(double const&, double const&, double const&)> Lqqb = [=] (double const& x, double const& z, double const& mu) -> double
    {
      double lqqb = 0;
      for (int i = 1; i <= apfel::NF(mu, Thresholds); i++)
	lqqb += apfel::QCh2[i-1] * ( ffs->xfxQ(i, z, mu) * pdfs->xfxQ(-i, x, mu) + ffs->xfxQ(-i, z, mu) * pdfs->xfxQ(i, x, mu) );
      return lqqb;
    };
  const std::function<double(double const&, double const&, double const&)> Lqqp1 = [=] (double const& x, double const& z, double const& mu) -> double
    {
      const int Nf = apfel::NF(mu, Thresholds);
      double lqqp1 = 0;
      for (int i = 1; i <= Nf; i++)
	for (int j = 1; j <= Nf; j++)
	  lqqp1 += apfel::QCh2[i-1] * ( ffs->xfxQ(j, z, mu) * pdfs->xfxQ(i, x, mu) + ffs->xfxQ(-j, z, mu) * pdfs->xfxQ(-i, x, mu) +
					ffs->xfxQ(j, z, mu) * pdfs->xfxQ(-i, x, mu) + ffs->xfxQ(-j, z, mu) * pdfs->xfxQ(i, x, mu) );
      return lqqp1;
    };
  const std::function<double(double const&, double const&, double const&)> Lqqp2 = [=] (double const& x, double const& z, double const& mu) -> double
    {
      const int Nf = apfel::NF(mu, Thresholds);
      double lqqp2 = 0;
      for (int i = 1; i <= Nf; i++)
	for (int j = 1; j <= Nf; j++)
	  lqqp2 += apfel::QCh2[j-1] * ( ffs->xfxQ(j, z, mu) * pdfs->xfxQ(i, x, mu) + ffs->xfxQ(-j, z, mu) * pdfs->xfxQ(-i, x, mu) +
					ffs->xfxQ(j, z, mu) * pdfs->xfxQ(-i, x, mu) + ffs->xfxQ(-j, z, mu) * pdfs->xfxQ(i, x, mu) );
      return lqqp2;
    };
  const std::function<double(double const&, double const&, double const&)> Lqqp3 = [=] (double const& x, double const& z, double const& mu) -> double
    {
      const int Nf = apfel::NF(mu, Thresholds);
      double lqqp3 = 0;
      for (int i = 1; i <= Nf; i++)
	for (int j = 1; j <= Nf; j++)
	  lqqp3 += apfel::QCh[i-1] * apfel::QCh[j-1] * ( ffs->xfxQ(j, z, mu) * pdfs->xfxQ(i, x, mu) + ffs->xfxQ(-j, z, mu) * pdfs->xfxQ(-i, x, mu) -
							 ffs->xfxQ(j, z, mu) * pdfs->xfxQ(-i, x, mu) - ffs->xfxQ(-j, z, mu) * pdfs->xfxQ(i, x, mu) );
      return lqqp3;
    };
  const std::function<double(double const&, double const&, double const&)> Lgq = [=] (double const& x, double const& z, double const& mu) -> double
    {
      double lgq = 0;
      for (int i = 1; i <= apfel::NF(mu, Thresholds); i++)
	lgq += apfel::QCh2[i-1] * ( pdfs->xfxQ(i, x, mu) + pdfs->xfxQ(-i, x, mu) );
      return ffs->xfxQ(21, z, mu) * lgq;
    };
  const std::function<double(double const&, double const&, double const&)> Lqg = [=] (double const& x, double const& z, double const& mu) -> double
    {
      double lqg = 0;
      for (int i = 1; i <= apfel::NF(mu, Thresholds); i++)
	lqg += apfel::QCh2[i-1] * ( ffs->xfxQ(i, z, mu) + ffs->xfxQ(-i, z, mu) );
      return lqg * pdfs->xfxQ(21, x, mu);
    };
  const std::function<double(double const&, double const&, double const&)> Lgg = [=] (double const& x, double const& z, double const& mu) -> double
    {
      return apfel::SumCh2[apfel::NF(mu, Thresholds)] * ffs->xfxQ(21, z, mu) * pdfs->xfxQ(21, x, mu);
    };

  // Define kinematics
  const double x = 0.01;
  const double z = 0.4;
  const double Q = 2;
  const int nf = apfel::NF(Q, Thresholds);

  // Now convolute parton luminosities with the coeffient functions
  // LO
  double FL = 0;
  double FT = LqqNS(x, z, Q);

  // NLO
  FL += pdfs->alphasQ(Q) / 2 / M_PI * (
    // qqNS
    + ( C1LQ2Q_Lx_Lz(nf) + C1LQ2Q_Lx_Sz(nf, z, true) + C1LQ2Q_Sx_Lz(nf, x, true) + C1LQ2Q_Sx_Sz(nf, x, z, true, true) ) * LqqNS(x, z, Q)                                                                                 // Lx_Lz
    + (apfel::Integrator  {[=] (                  double const& zp) -> double { return ( C1LQ2Q_Lx_Sz(nf, zp, false) + C1LQ2Q_Sx_Sz(nf, x, zp, true, false) ) * ( LqqNS(x, z / zp, Q) - LqqNS(x, z, Q) ); }}).integrate(z, 1, apfel::eps5) // Lx_Sz
    + (apfel::Integrator  {[=] (                  double const& zp) -> double { return ( C1LQ2Q_Lx_Rz(nf, zp) + C1LQ2Q_Sx_Rz(nf, x, zp, true) ) * LqqNS(x, z / zp, Q); }}).integrate(z, 1, apfel::eps5)                                    // Lx_Rz
    + (apfel::Integrator  {[=] (double const& xp                  ) -> double { return ( C1LQ2Q_Sx_Lz(nf, xp, false) + C1LQ2Q_Sx_Sz(nf, xp, z, false, true) ) * ( LqqNS(x / xp, z, Q) - LqqNS(x, z, Q) ); }}).integrate(z, 1, apfel::eps5) // Sx_Lz
    + (apfel::Integrator2D{[=] (double const& xp, double const& zp) -> double { return   C1LQ2Q_Sx_Sz(nf, xp, zp, false, false) * ( LqqNS(x / xp, z / zp, Q) - LqqNS(x, z / zp, Q) - LqqNS(x / xp, z, Q) + LqqNS(x, z, Q) ); }}).integrate(x, 1, z, 1, apfel::eps5) // Sx_Sz
    + (apfel::Integrator2D{[=] (double const& xp, double const& zp) -> double { return   C1LQ2Q_Sx_Rz(nf, xp, zp, false) * ( LqqNS(x / xp, z / zp, Q) - LqqNS(x, z / zp, Q) ); }}).integrate(x, 1, z, 1, apfel::eps5)      // Sx_Rz
    + (apfel::Integrator  {[=] (double const& xp                  ) -> double { return ( C1LQ2Q_Rx_Lz(nf, xp) + C1LQ2Q_Rx_Sz(nf, xp, z, true) ) * LqqNS(x / xp, z, Q); }}).integrate(z, 1, apfel::eps5)                                    // Rx_Lz
    + (apfel::Integrator2D{[=] (double const& xp, double const& zp) -> double { return   C1LQ2Q_Rx_Sz(nf, xp, zp, false) * ( LqqNS(x / xp, z / zp, Q) - LqqNS(x / xp, z, Q) ); }}).integrate(x, 1, z, 1, apfel::eps5)      // Rx_Sz
    + (apfel::Integrator2D{[=] (double const& xp, double const& zp) -> double { return   C1LQ2Q_Rx_Rz(nf, xp, zp) * LqqNS(x / xp, z / zp, Q); }}).integrate(x, 1, z, 1, apfel::eps5)                                       // Rx_Rz
    );

  // Print value
  std::cout << std::scientific << std::endl;
  std::cout << "Reduced SIDIS cross section (Q = " << Q << " GeV, x = " << x << ", z = " << z << "): " << CrossSectionFO(Q).Evaluate(x, z) << " - " << LqqNS(x, z, Q) << std::endl;

  return 0;
}
