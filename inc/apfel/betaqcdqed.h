//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name QCDxQED mixed beta function
   * Coefficients of the QCDxQED mixed \f$\beta\f$ function.
   * Reference: https://arxiv.org/pdf/hep-ph/9803211.pdf
   */
  ///@{
  /**
   * @brief \f$O(\alpha_s\alpha)\f$ coefficient of the QCD \f$\beta\f$
   * function.
   * @param nf: the number of active flavours
   * @return \f$\beta^{(\alpha_s\alpha)}(n_f)\f$
   */
  double beta1qcdqed(int const& nf);

  /**
   * @brief \f$O(\alpha\alpha_s)\f$ coefficient of the QCD \f$\beta\f$
   * function.
   * @param nf: the number of active flavours
   * @return \f$\beta^{(\alpha\alpha_s)}(n_f)\f$
   */
  double beta1qedqcd(int const& nf);
  ///@}
}
