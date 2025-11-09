/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Type aliases to model the nesting structure of OMEs
 */

#ifndef LIBOME_OME_TYPE_ALIASES_H
#define LIBOME_OME_TYPE_ALIASES_H

#include <apfel/ome/functions.h>
#include <apfel/ome/laurent_polynomial.h>
#include <apfel/ome/piecewise.h>

namespace apfel
{
  namespace ome
  {
    /**
     * \defgroup ome-data-types Public data types
     * \brief Specialised data types for the objects provided
     * @{
     */

    /**
     * \name Standard OMEs
     * \brief Type aliases for standard x-dependent OMEs, given by a piecewise
     *        definition based on generalised power series
     * @{
     */

    /**
     * \brief Truncated power series in \f$x-x_0\f$
     *
     * \details
     * Truncated power series that takes the value of \f$x-x_0\f$ and returns the
     * evaluation of the power series.
     */
    template<typename Tnum> using ome_x = laurent_polynomial<Tnum, Tnum>;

    /**
     * \brief Truncated power series in \f$\log(x-x_0)\f$ and \f$x-x_0\f$
     *
     * \details
     * Truncated power seires that takes the values of \f$x-x_0\f$ and
     * \f$\log(x-x_0)\f$ and returns the evaluation of the power series.
     */
    template<typename Tnum> using ome_logx
      = laurent_polynomial<Tnum, ome_x<Tnum>, Tnum>;

    /**
     * \brief Truncated generalised power series in x
     *
     * \details
     * Univariate wrapper for \ref ome_logx, that takes the value of \f$x\f$,
     * shifts it using the appropriate function, calculates \f$\log(x-x_0)\f$
     * and evaluates the underlying \ref ome_logx.
     */
    template<typename Tnum> using ome_genps
      = func_apply<Tnum, func_copy_and_log<Tnum, ome_logx<Tnum>>>;

    /**
     * \brief Construct an \ref ome_genps object
     *
     * \details
     * Since the template parameter deduction for \ref ome_genps using plain
     * constructors is tricky and require manual composition of \ref func_apply
     * and \ref func_copy_and_log classes, this helper allows to construct
     * \ref ome_genps objects with automatic template parameter deduction.
     *
     * \tparam Tnum Numerical type
     *
     * \param apply_function Univariate function that is applied to \f$x\f$ before
     *        it is passed to the wrapped \ref ome_logx target
     * \param target Underlying truncated power series in \f$\log(x-x_0)\f$ and
     *        \f$x-x_0\f$
     */
    template<typename Tnum>
    ome_genps<Tnum> make_ome_genps(std::function<Tnum(Tnum)> apply_function,
                                   ome_logx<Tnum> target)
    {
      return(
              func_apply<Tnum, func_copy_and_log<Tnum, ome_logx<Tnum>>>(
                apply_function,
                func_copy_and_log<Tnum, ome_logx<Tnum>>(target)
              )
            );
    };

    /**
     * \brief Piecewise defined OME coefficient in \f$x\f$
     *
     * \details
     * A collection of generalised power series representations along with the
     * intervals on which they are defined.
     */
    template<typename Tnum> using ome_piecewisex
      = piecewise<Tnum, ome_genps<Tnum>>;

    /**
     * \brief OME coefficients in \f$N_F\f$
     *
     * \details
     * Polynomial in the number of massless quarks \f$N_F\f$. Each coefficient is
     * a function of \f$x\f$.
     */
    template<typename Tnum> using ome_nf
      = laurent_polynomial<Tnum, ome_piecewisex<Tnum>, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$N_F\f$
     */
    template<typename Tnum> using ome_nf_view
      = laurent_polynomial_view<Tnum, ome_piecewisex<Tnum>, Tnum>;

    /**
     * \brief OME coefficients in \f$L_M\f$
     *
     * \details
     * Polynomial in the mass logarithm
     * \f$L_M = \log\left(\frac{m^2}{\mu^2}\right)\f$.
     * Each coefficient is a function of \f$N_F\f$ and \f$x\f$.
     */
    template<typename Tnum> using ome_logm
      = laurent_polynomial<Tnum, ome_nf<Tnum>, Tnum, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$L_M\f$
     */
    template<typename Tnum> using ome_logm_view
      = laurent_polynomial_view<Tnum, ome_nf<Tnum>, Tnum, Tnum>;

    /**
     * \brief OME coefficients in \f$a_s\f$
     *
     * \details
     * Polynomial in the strong coupling constant
     * \f$a_s(\mu) = \frac{\alpha_s(\mu)}{4 \pi} = \frac{g^2(\mu)}{(4 \pi)^2}\f$.
     * Each coefficient is a function of \f$L_M\f$, \f$N_F\f$ and \f$x\f$.
     *
     */
    template<typename Tnum> using ome_as
      = laurent_polynomial<Tnum, ome_logm<Tnum>, Tnum, Tnum, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$a_s\f$
     */
    template<typename Tnum> using ome_as_view
      = laurent_polynomial_view<Tnum, ome_logm<Tnum>, Tnum, Tnum, Tnum>;

    /**
     * @}
     */


    /**
     * \name Constant OMEs
     * \brief Type aliases for constant (x-independent) OMEs
     * @{
     */

    /**
     * \brief OME coefficients in \f$N_F\f$
     *
     * \details
     * Polynomial in the number of massless quarks \f$N_F\f$. The coefficients
     * are constant numbers.
     */
    template<typename Tnum> using ome_nf_const
      = laurent_polynomial<Tnum, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$N_F\f$
     */
    template<typename Tnum> using ome_nf_const_view
      = laurent_polynomial_view<Tnum, Tnum>;

    /**
     * \brief OME coefficients in \f$L_M\f$
     *
     * \details
     * Polynomial in the mass logarithm
     * \f$L_M = \log\left(\frac{m^2}{\mu^2}\right)\f$.
     * Each coefficient is a function of \f$N_F\f$.
     */
    template<typename Tnum> using ome_logm_const
      = laurent_polynomial<Tnum, ome_nf_const<Tnum>, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$L_M\f$
     */
    template<typename Tnum> using ome_logm_const_view
      = laurent_polynomial_view<Tnum, ome_nf_const<Tnum>, Tnum>;

    /**
     * \brief OME coefficients in \f$a_s\f$
     *
     * \details
     * Polynomial in the strong coupling constant
     * \f$a_s(\mu) = \frac{\alpha_s(\mu)}{4 \pi} = \frac{g^2(\mu)}{(4 \pi)^2}\f$.
     * Each coefficient is a function of \f$L_M\f$ and \f$N_F\f$.
     */
    template<typename Tnum> using ome_as_const
      = laurent_polynomial<Tnum, ome_logm_const<Tnum>, Tnum, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$a_s\f$
     */
    template<typename Tnum> using ome_as_const_view
      = laurent_polynomial_view<Tnum, ome_logm_const<Tnum>, Tnum, Tnum>;

    /**
     * @}
     */


    /**
     * \name Plus function part of OMEs
     * \brief Type aliases for plus function part of OMEs
     * @{
     */

    /**
     * \brief Coefficients of the basic plus functions
     *
     * \details
     * This polynomial contains the coefficients of the basic plus functions
     * \f[
     *   D_k(x) = \left[\frac{\log^k(1-x)}{1-x}\right]_+
     * \f]
     * In practice, it takes the numerical value of \f$\log(1-x)\f$ and
     * calculates the polynomial
     * \f[
     *   \sum_k \log^k(1-x) c_k
     * \f]
     * and the division by \f$(1-x)\f$ is done one level above in the hierarchy
     * in \ref ome_plusfunc. The coefficients \f$c_k\f$ are constant numbers.
     */
    template<typename Tnum> using ome_plusfunc_coeff
      = laurent_polynomial<Tnum, Tnum>;

    /**
     * \brief Sum of basic plus functions
     *
     * \details
     * This models a sum of basic plus distributions
     * \f[
     *   f_+(x) = \sum_k c_k \left[\frac{\log^k(1-x)}{1-x}\right]_+
     * \f]
     * It takes the argument x, shifts it to \f$(1-x)\f$, takes the logarithm
     * and passes it to the underlying \ref ome_plusfunc_coeff. The result is
     * then divided by \f$(1-x)\f$ before being returned.
     */
    template<typename Tnum> using ome_plusfunc
      = func_plusfunc_omx<Tnum, ome_plusfunc_coeff<Tnum>>;

    /**
     * \brief OME coefficients in \f$N_F\f$
     *
     * \details
     * Polynomial in the number of massless quarks \f$N_F\f$. The coefficients
     * are (plus) functions of \f$x\f$.
     */
    template<typename Tnum> using ome_nf_plus
      = laurent_polynomial<Tnum, ome_plusfunc<Tnum>, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$N_F\f$
     */
    template<typename Tnum> using ome_nf_plus_view
      = laurent_polynomial_view<Tnum, ome_plusfunc<Tnum>, Tnum>;

    /**
     * \brief OME coefficients in \f$L_M\f$
     *
     * \details
     * Polynomial in the mass logarithm
     * \f$L_M = \log\left(\frac{m^2}{\mu^2}\right)\f$.
     * Each coefficient is a function of \f$N_F\f$ and \f$x\f$.
     */
    template<typename Tnum> using ome_logm_plus
      = laurent_polynomial<Tnum, ome_nf_plus<Tnum>, Tnum, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$L_M\f$
     */
    template<typename Tnum> using ome_logm_plus_view
      = laurent_polynomial_view<Tnum, ome_nf_plus<Tnum>, Tnum, Tnum>;

    /**
     * \brief OME coefficients in \f$a_s\f$
     *
     * \details
     * Polynomial in the strong coupling constant
     * \f$a_s(\mu) = \frac{\alpha_s(\mu)}{4 \pi} = \frac{g^2(\mu)}{(4 \pi)^2}\f$.
     * Each coefficient is a function of \f$L_M\f$, \f$N_F\f$ and \f$x\f$.
     */
    template<typename Tnum> using ome_as_plus
      = laurent_polynomial<Tnum, ome_logm_plus<Tnum>, Tnum, Tnum, Tnum>;

    /**
     * \brief View type for OME coefficient in \f$a_s\f$
     */
    template<typename Tnum> using ome_as_plus_view
      = laurent_polynomial_view<Tnum, ome_logm_plus<Tnum>, Tnum, Tnum, Tnum>;

    /**
     * @}
     */

    /**
     * @}
     */
  }
}

#endif
