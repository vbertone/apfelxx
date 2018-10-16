//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Cusp anomalous dimension.
   * Coefficients of the &Gamma;<SUB>Cusp</SUB> anomalous
   * dimension. The expressions are taken from eq. (59)
   * https://arxiv.org/pdf/1705.07167.pdf.
   */
  ///@{
  /// &alpha;<SUB>s</SUB> term
  double GammaCusp0();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double GammaCusp1(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double GammaCusp2(int const& nf);
  ///@}
}
