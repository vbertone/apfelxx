/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

#ifndef LIBOME_POLAGGQ_H
#define LIBOME_POLAGGQ_H

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
     * \defgroup ome-data-polAggQ polAggQ
     * \brief Data objects for \f$\Delta A_{gg,Q}\f$
     * \{
     */
    /// Regular part of \f$\Delta A_{gg,Q}\f$
    extern const ome_as<double> polAggQ_reg;

    /// Plus part of \f$\Delta A_{gg,Q}\f$
    extern const ome_as_plus<double> polAggQ_plus;

    /// Delta part of \f$\Delta A_{gg,Q}\f$
    extern const ome_as_const<double> polAggQ_delta;

    /// Container for regular, plus and delta part of \f$\Delta A_{gg,Q}\f$
    extern const rpd_distribution<ome_as_view<double>, ome_as_plus_view<double>, ome_as_const_view<double>> polAggQ;

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
   * \defgroup c-interface-polAggQ polAggQ
   * \brief C interface functions for \f$\Delta A_{gg,Q}\f$
   * \{
   * \name Evaluation
   * \{
   */

/// Evaluate regular part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_reg(double as, double LM, double NF, double x);
/// Evaluate regular part of \f$\Delta A_{gg,Q}\f$ truncated in \f$a_s\f$
  double ome_polAggQ_reg_trunc_as(int trunc_order, double as, double LM, double NF, double x);
/// Evaluate a given \f$a_s\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_reg_coeff_as(int order_as, double LM, double NF, double x);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_reg_coeff_as_LM(int order_as, int order_LM, double NF, double x);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$, \f$N_F\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_reg_coeff_as_LM_NF(int order_as, int order_LM, int order_NF, double x);

/// Evaluate plus part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_plus(double as, double LM, double NF, double x);
/// Evaluate plus part of \f$\Delta A_{gg,Q}\f$ truncated in \f$a_s\f$
  double ome_polAggQ_plus_trunc_as(int trunc_order, double as, double LM, double NF, double x);
/// Evaluate a given \f$a_s\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_plus_coeff_as(int order_as, double LM, double NF, double x);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_plus_coeff_as_LM(int order_as, int order_LM, double NF, double x);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$, \f$N_F\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_plus_coeff_as_LM_NF(int order_as, int order_LM, int order_NF, double x);

/// Evaluate delta part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_delta(double as, double LM, double NF);
/// Evaluate delta part of \f$\Delta A_{gg,Q}\f$ truncated in \f$a_s\f$
  double ome_polAggQ_delta_trunc_as(int trunc_order, double as, double LM, double NF);
/// Evaluate a given \f$a_s\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_delta_coeff_as(int order_as, double LM, double NF);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_delta_coeff_as_LM(int order_as, int order_LM, double NF);
/// Evaluate a given \f$a_s\f$, \f$L_M\f$, \f$N_F\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  double ome_polAggQ_delta_coeff_as_LM_NF(int order_as, int order_LM, int order_NF);

  /**
   * \}
   *
   * \name Range getters
   * \{
   */

/// Get minimum power in \f$a_s\f$ of regular part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_reg_min_power();
/// Get maximum power in \f$a_s\f$ of regular part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_reg_max_power();
/// Get minimum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_reg_coeff_as_min_power(int order_as);
/// Get maximum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_reg_coeff_as_max_power(int order_as);
/// Get minimum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_reg_coeff_as_LM_min_power(int order_as, int order_LM);
/// Get maximum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the regular part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_reg_coeff_as_LM_max_power(int order_as, int order_LM);

/// Get minimum power in \f$a_s\f$ of plus part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_plus_min_power();
/// Get maximum power in \f$a_s\f$ of plus part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_plus_max_power();
/// Get minimum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_plus_coeff_as_min_power(int order_as);
/// Get maximum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_plus_coeff_as_max_power(int order_as);
/// Get minimum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_plus_coeff_as_LM_min_power(int order_as, int order_LM);
/// Get maximum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the plus part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_plus_coeff_as_LM_max_power(int order_as, int order_LM);

/// Get minimum power in \f$a_s\f$ of delta part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_delta_min_power();
/// Get maximum power in \f$a_s\f$ of delta part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_delta_max_power();
/// Get minimum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_delta_coeff_as_min_power(int order_as);
/// Get maximum power in \f$L_M\f$ of a specific \f$a_s\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_delta_coeff_as_max_power(int order_as);
/// Get minimum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_delta_coeff_as_LM_min_power(int order_as, int order_LM);
/// Get maximum power in \f$N_F\f$ of a specific \f$a_s\f$, \f$L_M\f$ coefficient of the delta part of \f$\Delta A_{gg,Q}\f$
  int ome_polAggQ_delta_coeff_as_LM_max_power(int order_as, int order_LM);

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
