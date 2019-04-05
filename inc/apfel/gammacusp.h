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
   * dimension. The expressions up to
   * O(&alpha;<SUB>s</SUB><SUP>3</SUP>) are taken from eq. (59)
   * https://arxiv.org/pdf/1705.07167.pdf. While the expressions at
   * are taken from O(&alpha;<SUB>s</SUB><SUP>4</SUP>) from eqs. (4.1)
   * and (4.2) of https://arxiv.org/pdf/1805.09638.pdf.
   * @note All the expressions do not include an overall factor
   * C<SUB>F</SUB> or C<SUB>A</SUB>.
   */
  ///@{
  /// &alpha;<SUB>s</SUB> term
  double GammaCusp0();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double GammaCusp1(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double GammaCusp2(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>4</SUP> term
  double GammaCusp3(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>4</SUP> correction to the gluon
  /// anonalous dimension (neglected for now).
  double GammaCusp3gmq(int const& nf);
  ///@}
}
