#include <iostream>
#include <iomanip>
#include <cmath>
#include <fstream>
#include "apfel/APFELevol.h"

int main()
{
  // Activate some options
  APFEL::SetFFNS(5);
  APFEL::SetPerturbativeOrder(1);
  //APFEL::SetAlphaEvolution("expanded");
  APFEL::InitializeAPFEL();

  const double nQ = 100;
  const double Qmin = 1;
  const double Qmax = 91.2;
  const double Qstp = exp( log( Qmax / Qmin ) / ( nQ - 1 ) );

  std::cout << std::scientific;

  APFEL::SetAlphaQCDRef(0.118, 91.2);
  for (double Q = Qmin; Q <= Qmax; Q *= Qstp)
    std::cout << Q << "\t" << APFEL::AlphaQCD(Q) << std::endl;
  std::cout << "\n";
  APFEL::SetAlphaQCDRef(APFEL::AlphaQCD(1), 1);
  for (double Q = Qmin; Q <= Qmax; Q *= Qstp)
    std::cout << Q << "\t" << APFEL::AlphaQCD(Q) << std::endl;

  return 0;
}
