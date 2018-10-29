//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name QED beta function
   * Coefficients of the QED \f$\beta_{QED}\f$ function.
   */
  ///@{
  /**
   * @brief LO coefficient of the QED \f$\beta\f$ function.
   * @param nf: the number of active quark flavours
   * @param nl: the number of active charged leptons
   * @return \f$\beta_0(n_f, n_l)\f$
   */
  double beta0qed(int const& nf, int const& nl);

  /**
   * @brief NLO coefficient of the QED \f$\beta\f$ function.
   * @param nf: the number of active flavours
   * @param nl: the number of active charged leptons
   * @return \f$\beta_1(n_f, n_l)\f$
   */
  double beta1qed(int const& nf, int const& nl);
  ///@}
}
