/*
  Author: Valerio Bertone
 */

// Standard libs
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <fstream>
#include <algorithm>
#include <glob.h>
#include <iomanip>
#include <locale>
#include <stdio.h>
#include <stdlib.h>
#include <functional>

// LHAPDF libs
#include "LHAPDF/LHAPDF.h"

// APFEL++ libs
#include <apfel/apfelxx.h>

using namespace std;
using namespace apfel;

// b* prescription
double bstar(double const& b, double const&)
{
  const double bmax = 2 * exp( - apfel::emc);
  return b / sqrt( 1 + pow(b / bmax, 2) );
}

// Non-perturnative function
double fNP(double const&, double const& b, double const& zetaf)
{
  const double g1 = 0.02;
  const double g2 = 0.5;
  const double Q0 = 1;
  return exp( ( - g1 - g2 * log( sqrt(zetaf) / Q0 / 2 ) ) * b * b / 2 );
}

double B2FOqq(double const& x, double const& z, double const& qT, double const& Q)
{
  return 2 * CF * ( ( 1 - x ) * ( 1 - z ) + 4 * x * z + ( 1 + x * x * z * z ) * Q * Q / x / z / qT / qT ); 
}

double B2FOqg(double const& x, double const& z, double const& qT, double const& Q)
{
  return 2 * TR * ( 8 * x * ( 1 - x ) + ( x * x + ( 1 - x ) * ( 1 - x ) ) * ( z * z + ( 1 - z ) * ( 1 - z ) ) * ( 1 - x ) * Q * Q / x / z / z / qT / qT ); 
}

double B2FOgq(double const& x, double const& z, double const& qT, double const& Q)
{
  return 2 * CF * ( ( 1 - x ) * z + 4 * x * ( 1 - z ) + ( 1 + x * x * ( 1 - z ) * ( 1 - z ) ) * ( 1 - z ) * Q * Q / x / z / z / qT / qT ); 
}

double BLFOqq(double const& x, double const& z)
{
  return 8 * CF * x * z;
}

double BLFOqg(double const& x, double const&)
{
  return 16 * TR * x * ( 1 -  x );
}

double BLFOgq(double const& x, double const& z)
{
  return 8 * CF * x * ( 1 - z );
}

int main() {
  // PDF and FF sets
  LHAPDF::PDF* PDFs = LHAPDF::mkPDF("CT14nlo");
  //LHAPDF::PDF* FFs  = LHAPDF::mkPDF("DSS07_NLO_HadronSum");
  LHAPDF::PDF* FFs  = LHAPDF::mkPDF("NNFF11_HadronSum_nlo");

  // Heavy quark masses.
  const double mc = PDFs->quarkThreshold(4);
  const double mb = PDFs->quarkThreshold(5);
  const double mt = PDFs->quarkThreshold(6);
  const vector<double> Masses{0, 0, 0, mc, mb, mt};

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
      const double Yp     = 1 + pow(1-y, 2);
      const int nf        = NF(Q, Masses);
      const double as     = PDFs->alphasQ(Q) / FourPi;
      const double alpha2 = pow(1./137,2);

      // Define integrands.
      const auto Xsec = [=] (double const& zb) -> double
      {
	const double xb0 = ( 1 - zb ) / ( 1 - zb * ( 1 - qToQ2 ) );
	const double xt  = x / xb0;
	const double zt  = z / zb;
	if (xt >= 1 || zt >=1)
	  return 0;

	// Distribution combinations
	double qq = 0;
	double qg = 0;
	double gq = 0;
	const double gPDF = PDFs->xfxQ(21,xt,Q);
	const double gFF  = FFs->xfxQ(21,zt,Q);
	for (int i = 1; i <= nf; i++)
	  {
	    qq += QCh2[i-1] * ( PDFs->xfxQ(i,xt,Q) * FFs->xfxQ(i,zt,Q) + PDFs->xfxQ(-i,xt,Q) * FFs->xfxQ(-i,zt,Q) );
	    qg += QCh2[i-1] * gPDF * ( FFs->xfxQ(i,zt,Q) + FFs->xfxQ(-i,zt,Q) );
	    gq += QCh2[i-1] * ( PDFs->xfxQ(i,xt,Q) + PDFs->xfxQ(-i,xt,Q) ) * gFF;
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
      const Integrator Integral{Xsec};
      return Integral.integrate(z, zmax, 1e-5);
    };

  // Initialize space- and time-like splitting functions.
  const Grid g{{SubGrid{50,1e-2,3}, SubGrid{60,1e-1,3}, SubGrid{50,5e-1,3}, SubGrid{60,7e-1,3}}};
  const auto PDFObj = InitializeDglapObjectsQCD(g, Masses);
  const auto FFObj  = InitializeDglapObjectsQCDT(g, Masses);

  // Ogata-quadrature object of degree zero.
  apfel::OgataQuadrature OgataObj2{0, 1e-11};
 
  // Function for the computation of the asymptotic cross section at
  // O(as).
  const auto xsecAsy = [=,&g] (double const& Q2, double const& x, double const& y, double const& z, double const& qT2) -> double
    {
      // Internal kinematic variables
      const double Q  = sqrt(Q2);
      const double qT = sqrt(qT2);

      const auto integrand1 = [=] (double const& b) -> double
	{
	  const double f1 = fNP(x, b, Q2);
	  const double f2 = f1;
	  const double lg = 2 * log(bstar(b, Q) * Q / 2 / exp( - apfel::emc));
	  return b * f1 * f2 * pow(lg, 1) / 2;
	};
      const auto integrand2 = [=] (double const& b) -> double
	{
	  const double f1 = fNP(x, b, Q2);
	  const double f2 = f1;
	  const double lg = 2 * log(bstar(b, Q) * Q / 2 / exp( - apfel::emc));
	  return b * f1 * f2 * pow(lg, 2) / 2;
	};

      // Useful definitions
      const double powsup = pow(qT/Q, 0.2);
      const double I1 = - ( 1 - powsup ) / qT2                     + powsup * OgataObj2.transform(integrand1, qT);
      const double I2 = - ( 1 - powsup ) * 2 * log(Q2 / qT2) / qT2 + powsup * OgataObj2.transform(integrand2, qT);
      const double S11 = 6 * CF;
      const double S12 = - 2 * CF;

      const double y2     = y * y;
      const double Yp     = 1 + pow(1-y, 2);
      const int nf        = NF(Q, Masses);
      const double as     = PDFs->alphasQ(Q) / FourPi;
      const double alpha2 = pow(1./137,2);

      // Create sets of distributions for PDFs and FFs.
      const map<int,Distribution> DistPDFs = DistributionMap(g, [=] (double const& x, double const& Q)->map<int,double>{ return PDFs->xfxQ(x,Q); }, Q);
      const map<int,Distribution> DistFFs  = DistributionMap(g, [=] (double const& z, double const& Q)->map<int,double>{ return FFs->xfxQ(z,Q); }, Q);

      DoubleObject<Distribution> DoubDist;
      for (int j = -nf; j <= nf; j++)
	{
	  if (j == 0)
	    continue;

	  const Set<Operator> Psl = PDFObj.at(nf).SplittingFunctions.at(0);
	  const Set<Operator> Ptl = FFObj.at(nf).SplittingFunctions.at(0);

	  const Distribution Ppdf = Psl.at(3) * DistPDFs.at(j) + Psl.at(4) * DistPDFs.at(21) / 2 / nf;
	  const Distribution Pff  = Ptl.at(3) * DistFFs.at(j)  + Ptl.at(4) * DistFFs.at(21) / 2 / nf;

	  DoubDist.AddTerm({QCh2[abs(j)-1] * ( S11 * I1 + S12 * I2 ), DistPDFs.at(j), DistFFs.at(j)});
	  DoubDist.AddTerm({ - QCh2[abs(j)-1] * I1, DistPDFs.at(j), Pff});
	  DoubDist.AddTerm({ - QCh2[abs(j)-1] * I1, Ppdf, DistFFs.at(j)});
	}

      return 2 * qT * 2 * M_PI * alpha2 * Yp * as * DoubDist.Evaluate(x,z) / z / x / y / Q2;
    };

  // Alpha_s (from PDFs)
  const auto Alphas = [&] (double const& mu) -> double{ return PDFs->alphasQ(mu); };

  // Rotate PDF and FF sets into the QCD evolution basis
  const auto RotPDFs = [=] (double const& x, double const& mu) -> std::map<int,double>{ return apfel::PhysToQCDEv(PDFs->xfxQ(x,mu)); };
  const auto RotFFs  = [=] (double const& x, double const& mu) -> std::map<int,double>{ return apfel::PhysToQCDEv(FFs->xfxQ(x,mu)); };

  // Construct set of distributions as a function of the scale to be
  // tabulated
  const auto EvolvedPDFs = [=,&g] (double const& mu) -> apfel::Set<apfel::Distribution>{
    return apfel::Set<apfel::Distribution>{apfel::EvolutionBasisQCD{apfel::NF(mu, Masses)}, DistributionMap(g, RotPDFs, mu)};
  };
  const auto EvolvedFFs = [=,&g] (double const& mu) -> apfel::Set<apfel::Distribution>{
    return apfel::Set<apfel::Distribution>{apfel::EvolutionBasisQCD{apfel::NF(mu, Masses)}, DistributionMap(g, RotFFs, mu)};
  };

  // Tabulate collinear PDFs and FFs
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabPDFs{EvolvedPDFs, 100, PDFs->qMin(), PDFs->qMax(), 3, Masses};
  const auto CollPDFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabPDFs.Evaluate(mu); };
  const apfel::TabulateObject<apfel::Set<apfel::Distribution>> TabFFs{EvolvedPDFs, 100, FFs->qMin(), FFs->qMax(), 3, Masses};
  const auto CollFFs = [&] (double const& mu) -> apfel::Set<apfel::Distribution> { return TabFFs.Evaluate(mu); };

  // Initialize TMD objects
  const auto TmdObj = apfel::InitializeTmdObjects(g, Masses);

  // Build evolved TMD PDFs and FFs
  const auto EvTMDPDFs = BuildTmdPDFs(TmdObj, CollPDFs, Alphas, 1, 1);
  const auto EvTMDFFs  = BuildTmdFFs(TmdObj, CollFFs, Alphas, 1, 1);

  // Ogata-quadrature object of degree zero (do not integrate in qT).
  apfel::OgataQuadrature OgataObj{};

  // Function for the computation of the asymptotic cross section at
  // O(as).
  const auto xsecRes = [=,&g] (double const& Q2, double const& x, double const& y, double const& z, double const& qT2) -> double
    {
      // Internal kinematic variables
      const double Q  = sqrt(Q2);
      const double qT = sqrt(qT2);

      // Useful definitions
      const double Yp     = 1 + pow(1 - y, 2);
      const int nf        = NF(Q, Masses);
      const double as     = PDFs->alphasQ(Q) / FourPi;
      const double alpha2 = pow(1./137,2);

      // Compute the hard factor
      const double hcs = apfel::HardFactorDY(1, as, nf, 1);

      // Construct the TMD luminosity in b space to be fed to be
      // trasformed in qT space.
      const auto TMDLumib = [=] (double const& b) -> double
      {
	// Get Evolved TMD PDFs and FFs and rotate them into the
	// physical basis
	const std::map<int,apfel::Distribution> xF = QCDEvToPhys(EvTMDPDFs(bstar(b, Q), Q, Q2).GetObjects());
	const std::map<int,apfel::Distribution> xD = QCDEvToPhys(EvTMDFFs(bstar(b, Q), Q, Q2).GetObjects());

	// Combine TMDs through the e.m. charges
	double lumi = 0;
	for (int i = 1; i <= nf; i++)
	  lumi += QCh2[i-1] * ( xF.at(i).Evaluate(x) * xD.at(i).Evaluate(z) + xF.at(-i).Evaluate(x) * xD.at(-i).Evaluate(z) );

	// Combine all pieces and return
	return b * lumi * fNP(x, b, Q2) * fNP(z, b, Q2);
      };

      return 2 * qT * 2 * M_PI * alpha2 * Yp * hcs * ( OgataObj.transform(TMDLumib, qT) / x / z ) / x / y / Q2;
    };

  // Kinematics in terms of the variables in which TIMBA is
  // differential.
  const double Q2 = 10;
  const double x  = 0.4;
  const double z  = 0.1;
  const double y  = 0.2;

  const int    nqT    = 20;
  const double qTmin  = 0.001;
  const double qTmax  = 3;
  const double qTstep = exp( log( qTmax / qTmin ) / ( nqT - 1 ) );

  cout << scientific;
  double qT = qTmin;
  for (int iqT = 0; iqT < nqT; iqT++)
    {
      cout << qT << "\t"
	   << xsecFO(Q2, x, y, z, qT) << "\t"
	   << xsecAsy(Q2, x, y, z, qT) << "\t"
	   << xsecRes(Q2, x, y, z, qT) << "\t"
	   << endl;
      qT *= qTstep;
    }

  return 0;
}
