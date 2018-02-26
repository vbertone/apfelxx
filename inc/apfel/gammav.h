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
   * Coefficients of the GammaV anomalous dimensions.
   */
  ///@{
  double gammaVq0();
  double gammaVq1(int const& nf);
  double gammaVq2(int const& nf);
  double gammaVg0(int const& nf);
  double gammaVg1(int const& nf);
  double gammaVg2(int const& nf);
  ///@}
}
