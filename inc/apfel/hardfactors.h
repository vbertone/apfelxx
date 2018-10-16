//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Hard factors
   */
  ///@{
 /**
   * @brief Perturbative hard factor for Drell-Yan.
   */
  double HardFactorDY(int const& pt, double const& Alphas, int const& nf, double const& kappa);

 /**
   * @brief Perturbative hard factor for SIDIS.
   */
  double HardFactorSIDIS(int const& pt, double const& Alphas, int const& nf, double const& kappa);
  ///@}
}
