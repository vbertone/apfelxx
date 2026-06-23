//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name GammaV anomalous dimension.  Coefficients of the
   * &gamma;<SUB>F</SUB> anomalous dimension. The expressions are
   * taken from eq. (58) of
   * https://arxiv.org/pdf/1705.07167.pdf. Coefficients of the
   * &gamma;<SUB>F</SUB> anomalous dimension at order
   * &alpha;<SUB>s</SUB><SUP>4</SUP> are taken from eqs. (16)-(17)
   * https://arxiv.org/pdf/2102.09725
   */
  ///@{
  /// Quark &alpha;<SUB>s</SUB> term
  double gammaFq0();

  /// Quark &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double gammaFq1(int const& nf);

  /// Quark &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double gammaFq2(int const& nf);

  /// Quark &alpha;<SUB>s</SUB><SUP>4</SUP> term
  double gammaFq3(int const& nf);

  /// Gluon &alpha;<SUB>s</SUB> term
  double gammaFg0(int const& nf);

  /// Gluon &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double gammaFg1(int const& nf);

  /// Gluon &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double gammaFg2(int const& nf);

  /// Gluon &alpha;<SUB>s</SUB><SUP>4</SUP> term
  double gammaFg3(int const& nf);
  ///@}
}
