#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "apfel/APFELevol.h"

int main()
{
  // Activate some options
  APFEL::SetPDFSet("MMHT2014nlo68cl");
  APFEL::SetAlphaQCDRef(0.12, 91.1876);
  APFEL::SetPoleMasses(1.4, 4.75, 172);
  APFEL::SetPerturbativeOrder(1);
  //APFEL::SetAlphaEvolution("expanded");
  //APFEL::SetPDFEvolution("truncated");
  APFEL::InitializeAPFEL();
  APFEL::EvolveAPFEL(100, 100);

  const double nx = 1000;
  const double xmin = 0.00001;
  const double xmax = 1;
  const double xstp = exp( log( xmax / xmin ) / ( nx - 1 ) );

  std::cout << std::scientific;
  for (double x = xmin; x <= xmax; x *= xstp)
    std::cout << x << "\t" << APFEL::xPDFj(0, x) << std::endl;

  return 0;
}
