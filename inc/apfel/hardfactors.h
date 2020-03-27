//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Hard factors
   * Coefficients of the hard function H for the available
   * processes.
   */
  ///@{
  /// &alpha;<SUB>s</SUB> correction to the Drell-Yan hard factor
  double H1DY();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> correction to the Drell-Yan hard factor
  double H2DY(int const& nf);

  /// &alpha;<SUB>s</SUB> correction to the SIDIS hard factor
  double H1SIDIS();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> correction to the SIDIS hard factor
  double H2SIDIS(int const& nf);

  /// &alpha;<SUB>s</SUB> correction to the H in gg fusion hard factor
  double H1ggH();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> correction to the H in gg fusion hard factor
  double H2ggH(int const& nf);
  ///@}
}
