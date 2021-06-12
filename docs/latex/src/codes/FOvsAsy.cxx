/*
  Author: Valerio Bertone
 */

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include <apfel/apfelxx.h>

double B2FOqq(double const& x, double const& z, double const& qT, double const& Q)
{
  return 2 * apfel::CF * ( ( 1 - x ) * ( 1 - z ) + 4 * x * z + ( 1 + x * x * z * z ) * Q * Q / x / z / qT / qT ); 
}

double B2FOqg(double const& x, double const& z, double const& qT, double const& Q)
{
  return 2 * apfel::TR * ( 8 * x * ( 1 - x ) + ( x * x + ( 1 - x ) * ( 1 - x ) ) * ( z * z + ( 1 - z ) * ( 1 - z ) ) * ( 1 - x ) * Q * Q / x / z / z / qT / qT ); 
}

double B2FOgq(double const& x, double const& z, double const& qT, double const& Q)
{
  return 2 * apfel::CF * ( ( 1 - x ) * z + 4 * x * ( 1 - z ) + ( 1 + x * x * ( 1 - z ) * ( 1 - z ) ) * ( 1 - z ) * Q * Q / x / z / z / qT / qT ); 
}

double BLFOqq(double const& x, double const& z)
{
  return 8 * apfel::CF * x * z;
}

double BLFOqg(double const& x, double const&)
{
  return 16 * apfel::TR * x * ( 1 -  x );
}

double BLFOgq(double const& x, double const& z)
{
  return 8 * apfel::CF * x * ( 1 - z );
}

int main() {
  // PDF and FF sets
  LHAPDF::PDF* PDFs = LHAPDF::mkPDF("NNPDF31_nlo_pch_as_0118");
  LHAPDF::PDF* FFs  = LHAPDF::mkPDF("MAPFF10NLOPIsum");

  // Heavy quark masses.
  const double mc = PDFs->quarkThreshold(4);
  const double mb = PDFs->quarkThreshold(5);
  const double mt = PDFs->quarkThreshold(6);
  const std::vector<double> Thresholds{0, 0, 0, mc, mb, mt};

  // Function for the computation of the fixed-order cross section at
  // O(as).
  const auto xsecFO = [=] (double const& Q2, double const& x, double const& y, double const& z, double const& qT2) -> double
    {
      // Internal kinematic variables
      const double Q  = sqrt(Q2);
      const double qT = sqrt(qT2);

      // Useful definitions
      const double qToQ2  = qT2 / Q2;
      const double zmax   = ( 1 - x ) / ( 1 - x * ( 1 - qToQ2 ) );
      const double y2     = y * y;
      const double Yp     = 1 + pow(1 - y, 2);
      const int nf        = apfel::NF(Q, Thresholds);
      const double as     = PDFs->alphasQ(Q) / apfel::FourPi;
      const double alpha2 = pow(1./137., 2);

      // Define integrands.
      const auto Xsec = [=] (double const& zb) -> double
      {
	const double xb0 = ( 1 - zb ) / ( 1 - zb * ( 1 - qToQ2 ) );
	const double xt  = x / xb0;
	const double zt  = z / zb;
	if (xt >= 1 || zt >= 1)
	  return 0;

	// apfel::Distribution combinations
	double qq = 0;
	double qg = 0;
	double gq = 0;
	const double gPDF = PDFs->xfxQ(21, xt, Q);
	const double gFF  = FFs->xfxQ(21, zt, Q);
	for (int i = 1; i <= nf; i++)
	  {
	    qq += apfel::QCh2[i-1] * ( PDFs->xfxQ(i,xt,Q) * FFs->xfxQ(i,zt,Q) + PDFs->xfxQ(-i,xt,Q) * FFs->xfxQ(-i,zt,Q) );
	    qg += apfel::QCh2[i-1] * gPDF * ( FFs->xfxQ(i,zt,Q) + FFs->xfxQ(-i,zt,Q) );
	    gq += apfel::QCh2[i-1] * ( PDFs->xfxQ(i,xt,Q) + PDFs->xfxQ(-i,xt,Q) ) * gFF;
	  }
	qq /= xt * zt;
	qg /= xt * zt;
	gq /= xt * zt;

	// Construct integrand
	const double integrand = as * x * xb0 *
	( ( B2FOqq(xb0, zb, qT, Q) - y2 / Yp * BLFOqq(xb0, zb) ) * qq
	  + ( B2FOqg(xb0, zb, qT, Q) - y2 / Yp * BLFOqg(xb0, zb) ) * qg
	  + ( B2FOgq(xb0, zb, qT, Q) - y2 / Yp * BLFOgq(xb0, zb) ) * gq ) / Q / Q / ( 1 - zb );

	return 2 * qT * 2 * M_PI * alpha2 * Yp * integrand / x / y / Q2;
      };

      // Integrate over z.
      const apfel::Integrator Integral{Xsec};
      return Integral.integrate(z, zmax, 1e-7);
    };

  // Initialize space- and time-like splitting functions.
  const apfel::Grid g{{{100, 1e-3, 3}, {60, 1e-1, 3}, {50, 5e-1, 3}, {60, 7e-1, 3}}};
  const auto PDFObj = InitializeDglapObjectsQCD(g, Thresholds);
  const auto FFObj  = InitializeDglapObjectsQCDT(g, Thresholds);
 
  // Function for the computation of the asymptotic cross section at
  // O(as).
  const auto xsecAsy = [=,&g] (double const& Q2, double const& x, double const& y, double const& z, double const& qT2) -> double
    {
      // Internal kinematic variables
      const double Q  = sqrt(Q2);
      const double qT = sqrt(qT2);

      // Useful definitions
      const double I1 = - 1 / qT2;
      const double I2 = - 2 * log(Q2 / qT2) / qT2;
      const double S11 = 6 * apfel::CF;
      const double S12 = - 2 * apfel::CF;

      const double y2     = y * y;
      const double Yp     = 1 + pow(1 - y, 2);
      const int nf        = apfel::NF(Q, Thresholds);
      const double as     = PDFs->alphasQ(Q) / apfel::FourPi;
      const double alpha2 = pow(1./137., 2);

      // Create sets of distributions for PDFs and FFs.
      const apfel::Set<apfel::Distribution> DistPDFs{apfel::EvolutionBasisQCD{apfel::NF(Q, Thresholds)}, apfel::DistributionMap(g, [=] (double const& x, double const& Q) -> std::map<int ,double>{ return apfel::PhysToQCDEv(PDFs->xfxQ(x, Q)); }, Q)};
      const apfel::Set<apfel::Distribution> DistFFs{apfel::EvolutionBasisQCD{apfel::NF(Q, Thresholds)}, apfel::DistributionMap(g, [=] (double const& x, double const& Q) -> std::map<int ,double>{ return apfel::PhysToQCDEv(FFs->xfxQ(x, Q)); }, Q)};

      const std::map<int, apfel::Distribution> xf = apfel::QCDEvToPhys(DistPDFs.GetObjects());
      const std::map<int, apfel::Distribution> zd = apfel::QCDEvToPhys(DistFFs.GetObjects());
      const std::map<int, apfel::Distribution> xPf = apfel::QCDEvToPhys((PDFObj.at(nf).SplittingFunctions.at(0) * DistPDFs).GetObjects());
      const std::map<int, apfel::Distribution> zPd = apfel::QCDEvToPhys((FFObj.at(nf).SplittingFunctions.at(0) * DistFFs).GetObjects());

      apfel::DoubleObject<apfel::Distribution> DoubDist;
      for (int j = -nf; j <= nf; j++)
	{
	  if (j == 0)
	    continue;

	  DoubDist.AddTerm({apfel::QCh2[abs(j)-1] * ( S11 * I1 + S12 * I2 ), xf.at(j), zd.at(j)});
	  DoubDist.AddTerm({ - apfel::QCh2[abs(j)-1] * I1, xf.at(j), zPd.at(j)});
	  DoubDist.AddTerm({ - apfel::QCh2[abs(j)-1] * I1, xPf.at(j), zd.at(j)});
	}

      return 2 * qT * 2 * M_PI * alpha2 * Yp * as * DoubDist.Evaluate(x,z) / z / x / y / Q2;
    };

  // Kinematics in terms of the variables in which TIMBA is
  // differential.
  const double EP  = 1;
  const double El  = 27.6;
  const double VS  = sqrt( 4 * EP * El );

  const double Q2 = 2;
  const double x  = 0.1;
  const double z  = 0.2;
  const double y  = Q2 / VS / VS / x;

  const int    nqT    = 1000;
  const double qTmin  = 0.0001;
  const double qTmax  = 1;
  const double qTstep = exp( log(qTmax / qTmin ) / ( nqT - 1 ) );

  std::cout << std::scientific;
  for (double qT = qTmin; qT <= qTmax; qT *= qTstep)
    std::cout << qT << "\t"
	      << xsecFO(Q2, x, y, z, qT) << "\t" << xsecAsy(Q2, x, y, z, qT) << "\t"
	      << std::endl;

  return 0;
}
