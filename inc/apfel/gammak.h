//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Cusp anomalous dimension.
   * Coefficients of the &gamma;<SUB>K</SUB> anomalous dimension. The
   * expressions up to O(&alpha;<SUB>s</SUB><SUP>3</SUP>) are taken
   * from eq. (59) https://arxiv.org/pdf/1705.07167.pdf. While the
   * expressions at O(&alpha;<SUB>s</SUB><SUP>4</SUP>) are taken from
   * Eq. (6.3) of https://arxiv.org/pdf/1911.10174.pdf or Eq. (6) of
   * https://arxiv.org/pdf/2002.04617v2.pdf
   * @note All the expressions do not include an overall factor
   * C<SUB>F</SUB> or C<SUB>A</SUB>.
   */
  ///@{
  /// &alpha;<SUB>s</SUB> term
  double gammaK0();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double gammaK1(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double gammaK2(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>4</SUP> term
  double gammaK3(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>4</SUP> correction to the quark
  /// anonalous dimension needed to obtain the gluon anomalous
  /// dimension (neglected for now).
  double gammaK3gmq(int const& nf);
  ///@}
}
