//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <cmath>
#include <vector>

namespace apfel
{
  /**
   * @defgroup NumericalConstants Numerical constants
   * Collection of numerical constants often used
   * in the computations.
   */
  ///@{
  ///@}
  /**
   * @defgroup MathConstants Mathematical constants
   * Collection of mathematical constants often used in the code.
   * @ingroup NumericalConstants
   */
  ///@{
  /**
   * @name Small numbers
   * @brief Small numbers used for cutoffs, integration accuracies, etc.
   */
  ///@{
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
  ///@}

  /**
   * @name Numerical constants
   * @brief Definitions for recurrent constants.
   */
  ///@{
  const double Pi2    = M_PI * M_PI;
  const double FourPi = 4 * M_PI;
  const double emc    = 0.5772156649015329;
  const double zeta2  = 1.6449340668482264; // Pi2 / 6;
  const double zeta3  = 1.2020569031595943;
  const double zeta4  = 1.0823232337111382; // Pi2 * Pi2 / 90;
  const double zeta5  = 1.0369277551433699;
  const double zeta6  = 1.0173430619844491;
  ///@}

  /**
   * @name QCD colour factors
   * The SU(3) Casimir's.
   */
  ///@{
  const double TR = 0.5;
  const double CF = 4. / 3.;
  const double CA = 3.;
  const double NC = 3.;
  ///@}
  ///@}

  /**
   * @defgroup PhysConstants Physical constants
   * Collection of physical constants often used in the code.
   * @ingroup NumericalConstants
   */
  ///@{
  /**
   * @name Quark charges
   * @brief Quark electric charges and their square.
   */
  ///@{
  const double ed  = - 1. / 3.;
  const double eu  =   2. / 3.;
  const double ed2 =   1. / 9.;
  const double eu2 =   4. / 9.;
  const std::vector<double> QCh  = {ed,  eu,  ed,  eu,  ed,  eu};
  const std::vector<double> QCh2 = {ed2, eu2, ed2, eu2, ed2, eu2};
  ///@}

  /**
   * @name Conversion factor
   * @brief Conversion factor from GeV<SUP>-2</SUP> to pb.
   */
  ///@{
  const double ConvFact = 0.3893793656e9;
  ///@}

  /**
   * @name Z-boson mass and width
   * @brief Value of the mass of the Z boson and its width in GeV
   * taken from:
   * http://pdg.lbl.gov/2018/listings/rpp2018-list-z-boson.pdf.
   */
  ///@{
  const double ZMass  = 91.1876;
  const double GammaZ = 2.4952;
  ///@}

  /**
   * @name W-boson mass and width
   * @brief Value of the mass of the W bosons and their width in GeV
   * taken from:
   * http://pdg.lbl.gov/2018/listings/rpp2018-list-w-boson.pdf.
   */
  ///@{
  const double WMass  = 80.379;
  const double GammaW = 2.085;
  ///@}

  /**
   * @name Proton mass
   * @brief Value of the mass of the proton in GeV taken from:
   * http://pdg.lbl.gov/2018/reviews/rpp2018-rev-phys-constants.pdf.
   */
  ///@{
  const double ProtonMass = 0.9382720813;
  ///@}

  /**
   * @name Weinberg angle
   * @brief Value of sin<SUP>2</SUP>&theta;<SUB>W</SUB> in the MSbar
   * scheme taken from:
   * http://pdg.lbl.gov/2018/reviews/rpp2018-rev-phys-constants.pdf.
   */
  ///@{
  const double Sin2ThetaW = 0.23122;
  ///@}

  /**
   * @name Fermi constant
   * @brief Value of G<SUB>F</SUB> in GeV<SUP>-2</SUP> taken from:
   * http://pdg.lbl.gov/2018/reviews/rpp2018-rev-phys-constants.pdf.
   */
  ///@{
  const double GFermi = 1.1663787e-5;
  ///@}

  /**
   * @name CKM matrix elements
   * @brief Absolute value of the CMK matrix elements and their square
   * taken from:
   * http://pdg.lbl.gov/2018/reviews/rpp2018-rev-ckm-matrix.pdf.
   */
  ///@{
  const double Vud  = 0.97446;
  const double Vus  = 0.22452;
  const double Vub  = 0.00365;
  const double Vcd  = 0.22438;
  const double Vcs  = 0.97359;
  const double Vcb  = 0.04214;
  const double Vtd  = 0.00896;
  const double Vts  = 0.04133;
  const double Vtb  = 0.999105;
  const double Vud2 = Vud * Vud;
  const double Vus2 = Vus * Vus;
  const double Vub2 = Vub * Vub;
  const double Vcd2 = Vcd * Vcd;
  const double Vcs2 = Vcs * Vcs;
  const double Vcb2 = Vcb * Vcb;
  const double Vtd2 = Vtd * Vtd;
  const double Vts2 = Vts * Vts;
  const double Vtb2 = Vtb * Vtb;
  const std::vector<double> CMK  = {Vud,  Vus,  Vub,  Vcd,  Vcs,  Vcb,  Vtd,  Vts,  Vtb};
  const std::vector<double> CKM2 = {Vud2, Vus2, Vub2, Vcd2, Vcs2, Vcb2, Vtd2, Vts2, Vtb2};
  ///@}
  ///@}
}
