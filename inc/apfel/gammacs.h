//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Collins-Soper anomalous dimension.
   * Coefficients of the Collins-Soper anomalous dimension. The
   * expressions are taken from eq. (69)
   * https://arxiv.org/pdf/1705.07167.pdf and from eq (D.9) of
   * https://arxiv.org/pdf/1604.07869.pdf.
   */
  ///@{
  /// &alpha;<SUB>s</SUB> term
  double CSd10();

/// &alpha;<SUB>s</SUB>L term
  double CSd11();

  /// &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double CSd20(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>2</SUP>L term
  double CSd21(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>2</SUP>L<SUP>2</SUP> term
  double CSd22(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double CSd30(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP>L term
  double CSd31(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP>L<SUP>2</SUP> term
  double CSd32(int const& nf);

  /// &alpha;<SUB>s</SUB><SUP>3</SUP>L<SUP>3</SUP> term
  double CSd33(int const& nf);
  ///@}
}
