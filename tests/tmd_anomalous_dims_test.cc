//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include <apfel/splittingfunctionsunp_sl.h>
#include <apfel/betaqcd.h>
#include <apfel/gammak.h>
#include <apfel/constants.h>

#include <iostream>

int main()
{
  const int nf = 5;

  const double A1qsf = apfel::P0ns{}.Singular(0);
  const double A2qsf = apfel::P1nsp{nf}.Singular(0);
  const double A3qsf = apfel::P2nsp{nf}.Singular(0);
  const double A4qsf = apfel::P3nsp{nf}.Singular(0);
  const double A1qad = apfel::CF * apfel::gammaK0() / 2;
  const double A2qad = apfel::CF * apfel::gammaK1(nf) / 2;
  const double A3qad = apfel::CF * apfel::gammaK2(nf) / 2;
  const double A4qad = apfel::CF * apfel::gammaK3(nf) / 2;

  std::cout << std::scientific;
  std::cout << "Ratio between quark cusp anomalous dimensions and coefficients of the singular part of the NS+ splitting functions:" << std::endl;
  std::cout << "LO   : " << A1qad / A1qsf << std::endl;
  std::cout << "NLO  : " << A2qad / A2qsf << std::endl;
  std::cout << "NNLO : " << A3qad / A3qsf << std::endl;
  std::cout << "NNNLO: " << A4qad / A4qsf << std::endl;
  std::cout << "\n";

  const double A1gsf = apfel::P0gg{nf}.Singular(0);
  const double A2gsf = apfel::P1gg{nf}.Singular(0);
  const double A3gsf = apfel::P2gg{nf}.Singular(0);
  const double A4gsf = apfel::P3gg{nf}.Singular(0);
  const double A1gad = apfel::CA * apfel::gammaK0() / 2;
  const double A2gad = apfel::CA * apfel::gammaK1(nf) / 2;
  const double A3gad = apfel::CA * apfel::gammaK2(nf) / 2;
  const double A4gad = apfel::CA * ( apfel::gammaK3(nf) + apfel::gammaK3gmq(nf) ) / 2;

  std::cout << std::scientific;
  std::cout << "Ratio between gluon cusp anomalous dimensions and coefficients of the singular part of the gluon-gluon splitting functions:" << std::endl;
  std::cout << "LO   : " << A1gad / A1gsf << std::endl;
  std::cout << "NLO  : " << A2gad / A2gsf << std::endl;
  std::cout << "NNLO : " << A3gad / A3gsf << std::endl;
  std::cout << "NNNLO: " << A4gad / A4gsf << std::endl;


  return 0;
}
