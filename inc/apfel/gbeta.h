//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name g-functions beta.
   * g-functions used for the analytic evolution of the strong coupling.
   */
  ///@{
  /**
   * @brief g-function for the LO analytic running of the strong
   * coupling.
   * @param lambda: evolution parameter
   * @return the g1 function
   */
  double g1beta(double const& lambda);

  /**
   * @brief g-function for the NLO analytic running of the strong
   * coupling.
   * @param nf: the number of active flavours
   * @param kappa: the resummation-scale parameter
   * @param lambda: evolution parameter
   * @return the g2 function
   */
  double g2beta(int const& nf, double const& kappa, double const& lambda);

  /**
   * @brief g-function for the NNLO analytic running of the
   * strong coupling.
   * @param nf: the number of active flavours
   * @param kappa: the resummation-scale parameter
   * @param lambda: evolution parameter
   * @return the g3 function
   */
  double g3beta(int const& nf, double const& kappa, double const& lambda);

  /**
   * @brief g-function for the NNNLO analytic running of the strong
   * coupling.
   * @param nf: the number of active flavours
   * @param kappa: the resummation-scale parameter
   * @param lambda: evolution parameter
   * @return the g3 function
   */
  double g4beta(int const& nf, double const& kappa, double const& lambda);
  ///@}
}
