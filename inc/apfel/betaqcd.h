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
   * @name QCD beta function
   * Coefficients of the QCD \f$\beta\f$ function.
   */
  ///@{
  /**
   * @brief LO coefficient of the QCD \f$\beta\f$ function.
   * @param nf: number of active flavours
   * @return \f$\beta_0(n_f)\f$
   */
  double beta0(int const& nf);

  /**
   * @brief NLO coefficient of the QCD \f$\beta\f$ function.
   * @param nf: number of active flavours
   * @return \f$\beta_1(n_f)\f$
   */
  double beta1(int const& nf);

  /**
   * @brief NNLO coefficient of the QCD \f$\beta\f$ function.
   * @param nf: number of active flavours
   * @return \f$\beta_2(n_f)\f$
   */
  double beta2(int const& nf);

  /**
   * @brief NNNLO coefficient of the QCD \f$\beta\f$ function.
   * @param nf: number of active flavours
   * @return \f$\beta_3(n_f)\f$
   */
  double beta3(int const& nf);
  ///@}
}
