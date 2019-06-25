/*
 * Author: Valerio Bertone
 */

// LHAPDF libs
#include <LHAPDF/LHAPDF.h>

// APFEL++ libs
#include <apfel/apfelxx.h>
#include <apfel/SIDIS.h>

int main() {
  // Open PDF and FF sets.
  const std::string PDFset = "CT14nlo";
  const std::string FFset  = "DSS14_NLO_PiSum";  // DSS14 set for pi^{\pm} converted in the LHAPDF format
  LHAPDF::PDF* PDFs = LHAPDF::mkPDF(PDFset);
  LHAPDF::PDF* FFs  = LHAPDF::mkPDF(FFset);

  // Define PDFs and FFs as lambda functions.
  auto fPDFs = [&] (double const& x, double const& Q) -> std::map<int, double> { return PDFs->xfxQ(x,Q); };
  auto fFFs  = [&] (double const& x, double const& Q) -> std::map<int, double> { return FFs->xfxQ(x,Q); };

  // Define x-space grid.
  const apfel::Grid g{{{100, 1e-4, 3}, {100, 1e-1, 3}, {50, 8e-1, 3}}};

  // Initialize SIDIS objects.
  InitializeSIDIS(g);

  // Define lambda function for the computation of the SIDIS
  // differential cross section dSigma/dxdQ2dz as a function of z, x,
  // Q2 and the squared center-of-mass energy of the collision S.
  const auto SigmaSIDIS = [&] (double const& x, double const& Q2, double const& z, double const& S) -> double
    {
      const double Q = sqrt(Q2);

      // Compute coupling from the LHAPDF PDF set.
      const double coup = PDFs->alphasQ(Q) / apfel::FourPi;

      // Compute distributions for PDFs and FFs.
      const std::map<int, apfel::Distribution> dPDF = DistributionMap(g, fPDFs, Q);
      const std::map<int, apfel::Distribution> dFF  = DistributionMap(g, fFFs,  Q);

      apfel::DoubleObject<apfel::Distribution> distqq;
      apfel::DoubleObject<apfel::Distribution> distgq;
      apfel::DoubleObject<apfel::Distribution> distqg;
      for (auto j = - 6; j <= 6; j++)
	{
	  // Skip the gluon.
	  if (j == 0 || dPDF.find(j) == dPDF.end() || dFF.find(j) == dFF.end())
	    continue;

	  distqq.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(j),  dFF.at(j)});
	  distgq.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(j),  dFF.at(21)});
	  distqg.AddTerm({apfel::QCh2[abs(j)-1], dPDF.at(21), dFF.at(j)});
	}

      // Inelasticity.
      const double y  = Q2 / S / x;
      const double y2 = - y * y;
      const double yp = 1 + pow(1 - y, 2);

      // DIS prefactor. ConvFact is the conversion factor from
      // GeV^{-2} to pb.
      const double fact = apfel::ConvFact * pow(1./137.,2) * apfel::FourPi / pow(Q, 4) / x;

      // Compute cross sections at LO and NLO.
      const apfel::DoubleObject<apfel::Distribution> xsec = yp * ( ( C20qq + coup * C21qq ) * distqq + coup * ( C21gq * distgq + C21qg * distqg ) )
      + y2 * coup * ( CL1qq * distqq + CL1gq * distgq + CL1qg * distqg );

      return fact * xsec.Evaluate(x, z);
    };

  // Kinematics
  const double x  = 0.157;
  const double Q2 = 20;
  const double z  = 0.35;
  const double y  = 0.439;
  const double S  = Q2 / x / y; // Q2 = xyS pow(318, 2); // HERA c.m.e.

  // Output cross section
  std::cout << std::scientific << "\ndSigma/dxdQ2dz(x = " << x << ", Q2 = " << Q2 << " GeV^2, z = " << z << ") = " << SigmaSIDIS(x, Q2, z, S) << " pb GeV^{-2}\n" << std::endl;

  return 0;
}
