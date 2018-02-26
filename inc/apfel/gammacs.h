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
   * @name Collins-Soper anomalous dimension.
   * Coefficients of the Collins-Soper anomalous dimension.
   */
  ///@{
  double CSd10();
  double CSd11();
  double CSd20(int const& nf);
  double CSd21(int const& nf);
  double CSd22(int const& nf);
  double CSd30(int const& nf);
  double CSd31(int const& nf);
  double CSd32(int const& nf);
  double CSd33(int const& nf);
  ///@}
}
