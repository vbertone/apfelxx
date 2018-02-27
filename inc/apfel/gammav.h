//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name GammaV anomalous dimension.
   * Coefficients of the &Gamma;<SUB>V</SUB> anomalous dimension.
   */
  ///@{
  /// Quark &alpha;<SUB>s</SUB> term
  double gammaVq0();
  /// Quark &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double gammaVq1(int const& nf);
  /// Quark &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double gammaVq2(int const& nf);
  /// Gluon &alpha;<SUB>s</SUB> term
  double gammaVg0(int const& nf);
  /// Gluon &alpha;<SUB>s</SUB><SUP>2</SUP> term
  double gammaVg1(int const& nf);
  /// Gluon &alpha;<SUB>s</SUB><SUP>3</SUP> term
  double gammaVg2(int const& nf);
  ///@}
}
