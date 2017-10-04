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
  const double emc    = 0.5772156649;
  const double zeta2  = 1.644934067;
  const double zeta3  = 1.2020569031;
  const double zeta4  = 1.0823232337;
  const double zeta5  = 1.0369277551;

  // Quark electric charges and their square
  const double ed  = - 1. / 3.;
  const double eu  =   2. / 3.;
  const double ed2 =   1. / 9.;
  const double eu2 =   4. / 9.;
  const vector<double> QCh  = {ed, eu, ed, eu, ed, eu};
  const vector<double> QCh2 = {ed2, eu2, ed2, eu2, ed2, eu2};

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

  // Zero's of the Bessel function j0.
  const vector<double> j0Zeros{2.4048255576957730, 5.5200781102863115, 8.6537279129110130,
      11.791534439014281, 14.930917708487787, 18.071063967910924,
      21.211636629879260, 24.352471530749302, 27.493479132040253,
      30.634606468431976, 33.775820213573570, 36.917098353664045,
      40.058425764628240, 43.199791713176730, 46.341188371661815,
      49.482609897397815, 52.624051841115000, 55.765510755019980,
      58.906983926080940, 62.048469190227166, 65.189964800206870,
      68.331469329856800, 71.472981603593740, 74.614500643701830,
      77.756025630388050, 80.897555871137630, 84.039090776938200,
      87.180629843641160, 90.322172637210490, 93.463718781944780,
      96.605267950996260, 99.746819858680600, 102.88837425419480,
      106.02993091645162, 109.17148964980538, 112.31305028049490,
      115.45461265366694, 118.59617663087253, 121.73774208795096,
      124.87930891323295, 128.02087700600833, 131.16244627521390,
      134.30401663830546, 137.44558802028428, 140.58716035285428,
      143.72873357368974, 146.87030762579664, 150.01188245695477,
      153.15345801922788, 156.29503426853353, 159.43661116426316,
      162.57818866894667, 165.71976674795502, 168.86134536923583,
      172.00292450307820, 175.14450412190274, 178.28608420007376,
      181.42766471373105, 184.56924564063870, 187.71082696004936,
      190.85240865258152, 193.99399070010912, 197.13557308566140,
      200.27715579333240};

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
