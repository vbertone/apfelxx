/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#ifndef LIBOME_AGQQ_H
#define LIBOME_AGQQ_H

#ifdef __cplusplus

#include <apfel/ome/ome_type_aliases.h>
#include <apfel/ome/rpd_distribution.h>

namespace apfel
{
  namespace ome
  {
    /**
     * \addtogroup ome-data
     * \{
     * \defgroup ome-data-AgqQ AgqQ
     * \brief Data objects for \f$A_{gq,Q}\f$
     * \{
     */
    /// Regular part of \f$A_{gq,Q}\f$
    extern const ome_as<double> AgqQ_reg;

    /// Container for regular, plus and delta part of \f$A_{gq,Q}\f$
    extern const rpd_distribution<ome_as_view<double>, ome_as_plus_view<double>, ome_as_const_view<double>> AgqQ;

    /**
     * \}
     * \}
     */
  }

  extern "C"
  {
#endif /* ifdef __cpluscplus */

  /* C interface */
  /**
   * \addtogroup c-interface
   * \{
   * \defgroup c-interface-AgqQ AgqQ
   * \brief C interface functions for \f$A_{gq,Q}\f$
   * \{
   * \name Evaluation
   * \{
   */

/// Evaluate regular part of \f$A_{gq,Q}\f$
  double ome_AgqQ_reg(double as, double LM, double NF, double x);
/// Evaluate regular part of \f$A_{gq,Q}\f$ truncated in \f$a_s\f$
  double ome_AgqQ_reg_trunc_as(int trunc_order, double as, double LM, double NF, double x);
/// Evaluate a given \f$a_s\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  double ome_AgqQ_reg_coeff_as(int order_as, double LM, double NF, double x);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  double ome_AgqQ_reg_coeff_as_LM(int order_as, int order_LM, double NF, double x);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$, \f$N_F\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  double ome_AgqQ_reg_coeff_as_LM_NF(int order_as, int order_LM, int order_NF, double x);



  /**
   * \}
   *
   * \name Range getters
   * \{
   */

/// Get minimum power in \f$a_s\f$ of regular part of \f$A_{gq,Q}\f$
  int ome_AgqQ_reg_min_power();
/// Get maximum power in \f$a_s\f$ of regular part of \f$A_{gq,Q}\f$
  int ome_AgqQ_reg_max_power();
/// Get minimum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  int ome_AgqQ_reg_coeff_as_min_power(int order_as);
/// Get maximum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  int ome_AgqQ_reg_coeff_as_max_power(int order_as);
/// Get minimum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  int ome_AgqQ_reg_coeff_as_LM_min_power(int order_as, int order_LM);
/// Get maximum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the regular part of \f$A_{gq,Q}\f$
  int ome_AgqQ_reg_coeff_as_LM_max_power(int order_as, int order_LM);



  /**
   * \}
   * \}
   * \}
   */

#ifdef __cplusplus
}
}
#endif
#endif
