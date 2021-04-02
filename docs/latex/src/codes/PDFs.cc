#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "apfel/APFELevol.h"
#include "LHAPDF/LHAPDF.h"

extern "C" void externalsetapfel_(double const& x, double const& Q, double* xf);

int main()
{
  // Activate some options
/*
  // MMHT
  APFEL::SetPDFSet("MMHT2014nlo68cl");
  APFEL::SetAlphaQCDRef(0.12, 91.1876);
  APFEL::SetPoleMasses(1.4, 4.75, 172);
  APFEL::SetPerturbativeOrder(1);
  //APFEL::SetAlphaEvolution("expanded");
  //APFEL::SetPDFEvolution("truncated");
  APFEL::InitializeAPFEL();
  APFEL::EvolveAPFEL(100, 100);

  // NNPDF
  APFEL::SetPDFSet("NNPDF31_nlo_pch_as_0118");
  APFEL::SetAlphaQCDRef(0.118, 91.2);
  APFEL::SetPoleMasses(1.51, 4.92, 172);
  APFEL::SetPerturbativeOrder(1);
  //APFEL::SetAlphaEvolution("expanded");
  //APFEL::SetPDFEvolution("truncated");
  APFEL::InitializeAPFEL();
  APFEL::EvolveAPFEL(100, 100);
*/
  const double nx = 1000;
  const double xmin = 0.00001;
  const double xmax = 1;
  const double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );
/*
  std::cout << std::scientific;
  for (double x = xmin; x <= xmax; x *= xstp)
    std::cout << x << "\t" << APFEL::xPDFj(0, x) << std::endl;
  APFEL::CleanUp();
*/
  // Hysteresis with MMHT
  const double Q1 = 5;
  const double Q2 = 100;
  APFEL::SetAlphaQCDRef(0.12, 91.1876);
  APFEL::SetPoleMasses(1.4, 4.75, 172);
  APFEL::SetPerturbativeOrder(1);
  //APFEL::SetAlphaEvolution("expanded");
  //APFEL::SetPDFEvolution("truncated");
  APFEL::InitializeAPFEL();
  APFEL::SetPDFSet("MMHT2014nlo68cl");

  APFEL::EvolveAPFEL(Q1, Q2);
  APFEL::SetPDFSet("external");

  APFEL::EvolveAPFEL(Q2, Q1);
  const LHAPDF::PDF* lhpdf = LHAPDF::mkPDF("MMHT2014nlo68cl");

  std::cout << std::scientific;
  for (double x = xmin; x <= xmax; x *= xstp)
    std::cout << x << "\t" << APFEL::xPDFj(0, x) << "\t" << lhpdf->xfxQ(21, x, Q1) << std::endl;

  return 0;
}

void externalsetapfel_(double const& x, double const& Q, double* xf)
{
  for (int i = -6; i <= 6; i++)
    xf[i+6] = APFEL::xPDFj(i, x);
}
