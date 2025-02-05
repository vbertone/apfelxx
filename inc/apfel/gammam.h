//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name MSbar mass anomalous dimension
   * Coefficients of the \f$\gamma_m\f$ function.
   * Reference: https://arxiv.org/pdf/hep-ph/9910332
   */
  ///@{
  /**
   * @brief LO coefficient of the \f$\gamma_m\f$ function.
   * @return \f$\gamma_m^{[0]}\f$
   */
  double gammam0();

  /**
   * @brief NLO coefficient of the \f$\\gamma_m\f$ function.
   * @param nf: the number of active flavours
   * @return \f$\gamma_m^{[1]}(n_f)\f$
   */
  double gammam1(int const& nf);

  /**
   * @brief NNLO coefficient of the \f$\\gamma_m\f$ function.
   * @param nf: the number of active flavours
   * @return \f$\gamma_m^{[2]}(n_f)\f$
   */
  double gammam2(int const& nf);
  ///@}
}
