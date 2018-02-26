//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <cmath>
#include <vector>

using namespace std;

namespace apfel
{
  /**
   * @brief Small numbers used for cutoffs, integration accuracies,
   * and so on.
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

  /**
   * @brief QCD colour factor (i.e. SU(3)).
   */
  const double TR = 0.5;
  const double CF = 4./3.;
  const double CA = 3.;

  /**
   * @brief Numerical constants.
   */
  const double Pi2    = M_PI * M_PI;
  const double FourPi = 4 * M_PI;
  const double emc    = 0.5772156649015329;
  const double zeta2  = 1.6449340668482264; // Pi2 / 6;
  const double zeta3  = 1.2020569031595943;
  const double zeta4  = 1.0823232337111382; // Pi2 * Pi2 / 90;
  const double zeta5  = 1.0369277551433699;

  /**
   * @brief Quark electric charges and their square.
   */
  const double ed  = - 1. / 3.;
  const double eu  =   2. / 3.;
  const double ed2 =   1. / 9.;
  const double eu2 =   4. / 9.;
  const vector<double> QCh  = {ed, eu, ed, eu, ed, eu};
  const vector<double> QCh2 = {ed2, eu2, ed2, eu2, ed2, eu2};

  /**
   * @brief Conversion factor from GeV^{-2} to pb.
   */
  const double ConvFact = 0.389379e9;
}
