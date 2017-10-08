//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <iostream>
#include <exception>
#include <stdexcept>
#include <string>
#include <sstream>
#include <vector>
#include <cmath>

using namespace std;

namespace apfel
{
  /**
    * Some general tools for the output and exception traitment.
    */
  const double eps2  = 1e-2;
  const double eps3  = 1e-3;
  const double eps4  = 1e-4;
  const double eps5  = 1e-5;
  const double eps6  = 1e-6;
  const double eps7  = 1e-7;
  const double eps8  = 1e-8;
  const double eps9  = 1e-9;
  const double eps10 = 1e-10;
  const double eps11 = 1e-11;
  const double eps12 = 1e-12;
  const double eps13 = 1e-13;
  const double eps14 = 1e-14;
  const double eps15 = 1e-15;
  const double eps25 = 1e-25;

  // Color factors
  const double TR = 0.5;
  const double CF = 4./3.;
  const double CA = 3.;

  // Useful constants
  const double Pi2    = M_PI * M_PI;
  const double FourPi = 4 * M_PI;
  const double emc    = 0.5772156649015329;
  const double zeta2  = 1.6449340668482264; // Pi2 / 6;
  const double zeta3  = 1.2020569031595943;
  const double zeta4  = 1.0823232337111382; // Pi2 * Pi2 / 90;
  const double zeta5  = 1.0369277551433699;

  // Quark electric charges and their square
  const double ed  = - 1. / 3.;
  const double eu  =   2. / 3.;
  const double ed2 =   1. / 9.;
  const double eu2 =   4. / 9.;
  const vector<double> QCh  = {ed, eu, ed, eu, ed, eu};
  const vector<double> QCh2 = {ed2, eu2, ed2, eu2, ed2, eu2};

  // Conversion factor from GeV^{-2} to pb
  const double ConvFact = 0.389379e9;

  // QCD Beta functions coefficients
  double beta0(int const& nf);
  double beta1(int const& nf);
  double beta2(int const& nf);

  // Cusp anomalous dimension coefficients
  double GammaCusp0();
  double GammaCusp1(int const& nf);
  double GammaCusp2(int const& nf);

  // Anomalous dimension gammaV coefficients
  double gammaVq0();
  double gammaVq1(int const& nf);
  double gammaVq2(int const& nf);
  double gammaVg0(int const& nf);
  double gammaVg1(int const& nf);
  double gammaVg2(int const& nf);

  // Collins-Soper anomalous dimension coefficients
  double CSd10();
  double CSd11();
  double CSd20(int const& nf);
  double CSd21(int const& nf);
  double CSd22(int const& nf);
  double CSd30(int const& nf);
  double CSd31(int const& nf);
  double CSd32(int const& nf);
  double CSd33(int const& nf);

  // Coeffiecient of log(zeta) in the zeta-prescription scheme.
  double Lzetaq10();
  double Lzetaq11();
  double Lzetaq20(int const& nf);
  double Lzetaq22(int const& nf);
  double Lzetag10(int const& nf);
  double Lzetag11();
  double Lzetag20(int const& nf);
  double Lzetag22(int const& nf);

  /**
   * @brief Return the number of active flavours at the scale Q given
   * the (ordered) vector of thresholds.
   * @param Q scale
   * @param Thresholds
   * @return number of active flavours
   */
  int NF(double const& Q, vector<double> const& Thresholds);

  //! Enumerator for code warning
  enum code {red = 31, green = 32, yellow = 33, blue = 34, normal = 39};

  /**
   * @brief info
   * @param tag
   * @param what
   */
  void info(string const& tag, string const& what);

  /**
   * @brief warning
   * @param tag
   * @param what
   */
  void warning(string const& tag, string const& what);

  /**
   * @brief success
   * @param tag
   * @param what
   */
  void success(string const& tag, string const& what);

  /**
   * @brief error
   * @param tag
   * @param what
   * @return
   */
  string error(string const& tag, string const& what);

  /**
   * @brief The runtime_exception class
   */
  class runtime_exception: public runtime_error
  {
  public:
    runtime_exception(string const& tag,
                      string const& what):
      runtime_error(error(tag,what)) {}
  };

  /**
   * @brief The logic_exception class
   */
  class logic_exception: public logic_error
  {
  public:
    logic_exception(string const& tag,
                    string const& what):
      logic_error(error(tag,what)) {}
  };
}

namespace std {
  /**
   * @brief implementation of operator<< for apfel::code
   */
  ostream& operator<<(ostream& os, apfel::code code);
}
