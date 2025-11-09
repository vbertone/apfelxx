/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Function wrappers that modify arguments passed to a wrapped callable
 */

#ifndef LIBOME_FUNCTIONS_H
#define LIBOME_FUNCTIONS_H

#include <type_traits>
#include <functional>
#include <algorithm>
#include <numeric>
#include <limits>
#include <cmath>
#include <cassert>

namespace apfel
{
  namespace ome
  {
    /// Single parameter identity function
    template<typename Tnum>
    Tnum f_id(Tnum x) { return(x); }

    /// Single parameter function: 1/2-x
    template<typename Tnum>
    Tnum f_half(Tnum x) { return(static_cast<Tnum>(1)/static_cast<Tnum>(2) - x); }

    /// Single parameter function: 1-x
    template<typename Tnum>
    Tnum f_omx(Tnum x) { return(static_cast<Tnum>(1)-x); }

    /// Single parameter function: (a+b*x)
    template<typename Tnum, int a, int b>
    Tnum f_linear(Tnum x)
    {
      return(static_cast<Tnum>(a) + static_cast<Tnum>(b)*x);
    }

    /// Single parameter function: (a+b*x)/(c+d*x)
    template<typename Tnum, int a, int b, int c, int d>
    Tnum f_moebius(Tnum x)
    {
      return((static_cast<Tnum>(a) + static_cast<Tnum>(b)*x) /
             (static_cast<Tnum>(c) + static_cast<Tnum>(d)*x));
    }

    /**
     * \brief Wrapper for callable that shifts the first argument by a fixed value
     *
     * \details
     * Takes a shift and callable at construction and acts as a wrapper which,
     * upon evaluation, shifts the first argument by the shift.
     *
     * \tparam Tnum Numerical type for the shift, the first argument and return
     *         type of the evaluation operator.
     * \tparam Tfunc Callable type to wrap
     * \tparam Trest Type parameter pack for the remainder of the arguments
     */
    template<typename Tnum, typename Tfunc, typename... Trest>
    class func_shift
    {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the callable type template parameter
      using element_type = Tfunc;

      /**
       * \brief Default constructor
       *
       * \details
       * Calls the default constructor on the target function class and
       * initialises shift to zero (i.e. the func_shift class acts as an
       * identity).
       */
      func_shift()
        : target_function_(element_type()),
          shift_(static_cast<numeric_type>(0)) {};

      /**
       * \brief Construct wrapper with given shift and callable
       *
       * \param shift Value of the shift
       * \param target_function Callable object to wrap
       */
      func_shift(numeric_type shift, element_type target_function)
        : target_function_(target_function), shift_(shift) {};

      /**
       * \brief Evaluate the wrapped callable with shifted first argument
       *
       * \param x Argument that gets shifted
       * \param rest Function parameter pack that is passed on to the wrapped
       *        callable unchanged
       *
       * \return The value returned by the wrapped callable
       */
      numeric_type operator()(numeric_type x, Trest... rest) const
      {
        return(target_function_(x+shift_, rest...));
      };

    private:
      element_type target_function_;
      numeric_type shift_;
    };

    /**
     * \brief Wrapper for callable that applies a function to the first argument
     *
     * \details
     * Wraps a callable and applies a function to the first argument upon
     * evaluation. The callalbe and function are specified at the time of
     * construction.
     *
     * \tparam Tnum Numerical type for the first argument and return type
     *         of the evaluation operator
     * \tparam Tfunc Callable type to wrap
     * \tparam Trest Pack of types for the remaining arguments of the wrapped
     *         callable
     */
    template<typename Tnum, typename Tfunc, typename... Trest>
    class func_apply
    {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the callable type template parameter
      using element_type = Tfunc;

      /**
       * \brief Default constructor
       *
       * \details
       * Calls default constructor on the target function and initialises the
       * apply function to f_id<numeric_type> (that is func_apply behaves as
       * the identity).
       */
      func_apply()
        : target_function_(element_type()),
          apply_function_(f_id<numeric_type>) {};

      /**
       * \brief Construct wrapper with given function and callable
       *
       * \param apply_function Function to apply to the first agument
       * \param target_function Callable object to wrap
       */
      func_apply(std::function<numeric_type(numeric_type)> apply_function,
                 element_type target_function)
        : target_function_(target_function),
          apply_function_(apply_function) {};

      /**
       * \brief Evaluate the wrapped callable with the first argument fed
       *        through the specified function
       *
       * \param x Argument that gets modified by the function
       * \param rest Function parameter pack that is passed on to the wrapped
       *        callable unchanged
       *
       * \return The value returned by the wrapped callable
       */
      numeric_type operator()(numeric_type x, Trest... rest) const
      {
        return(target_function_(apply_function_(x), rest...));
      };

    private:
      element_type target_function_;
      std::function<numeric_type(numeric_type)> apply_function_;
    };

    /**
     * \brief Wrapper for callable that copies the first argument and applies
     *        the logarithm to the first copy
     *
     * \details
     * Wraps a callable that takes one more argument than this wrapper. Upon
     * evaluation this wrapper passes the first argument of the wrapper to
     * the first two arguments of the wrapped callable. Moreover, it applies
     * the logarithm to the first passed argument. I.e. if f wraps g and f is
     * called as f(x,...) it calls g(log(x),x,...).
     *
     * \tparam Tnum Numerical type for the first argument and the return type of
     *         the wrapped callable
     * \tparam Tfunc Callable type to wrap
     * \tparam Trest Pack of types for the remaining arguments of the wrapped
     *         callable
     */
    template<typename Tnum, typename Tfunc, typename... Trest>
    class func_copy_and_log
    {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the callable type template parameter
      using element_type = Tfunc;

      /**
       * \brief Default constructor
       *
       * \details
       * Calls the default constructor on the target function.
       */
      func_copy_and_log()
        : target_function_(element_type()) {};

      /**
       * \brief Construct wrapper from callable
       *
       * \param target_function Callable object to wrap
       */
      explicit
      func_copy_and_log(element_type target_function)
        : target_function_(target_function) {};

      /**
       * \brief Evaluate the wrapped callable object with the first argument
       *        copied and with log() applied to the first copy
       *
       * \param x Argument to copy
       * \param rest Function parameter pack that is passed on to the wrapped
       *        callable unchanged
       *
       * \return The value returned by the wrapped callable
       */
      numeric_type operator()(numeric_type x, Trest... rest) const
      {
        using std::log;
        return(target_function_(x > static_cast<numeric_type>(0) ? log(x)
                                : std::numeric_limits<numeric_type>::quiet_NaN(), x, rest...));
      };

    private:
      element_type target_function_;
    };

    /**
     * \brief Wrapper for callable that emulates \f$1-x\f$ plus function kernels
     *
     * \details
     * Upon evaluation this wrapper takes the first argument \f$x\f$, takes the
     * logarithm \f$\log(1-x)\f$ of it and passes it on to the wrapped callable.
     * The result from the callable is then divided by \f$1-x\f$. If \f$f\f$
     * wraps \f$g\f$, evaluating \f$f\f$ corresponds to
     * \f[
     *   f(x,\dots) = \frac{g(\log(1-x),\dots)}{1-x}
     * \f]
     * The idea is that by wrapping a laurent_polynomial with this wrapper, it is
     * straightforward to construct a sum of plus function kernels
     * \f[
     *   p(x,\dots) = \sum_k \frac{\log^k(1-x)}{1-x} c_i(\dots)
     * \f]
     *
     * \tparam Tnum Numerical type for the first argument and the return type of
     *         the callable
     * \tparam Tfunc Callable type to wrap
     * \tparam Trest Pack of types for the remaining arguments of the wrapped
     *         callable
     */
    template<typename Tnum, typename Tfunc, typename... Trest>
    class func_plusfunc_omx
    {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the callable type template parameter
      using element_type = Tfunc;
      /// Boolean type alias indicating that this class has an eval_plus_int method
      using has_eval_plus_int = std::true_type;

      /**
       * \brief Default constructor
       *
       * \details
       * Calls the default constructor on the target function class.
       */
      func_plusfunc_omx()
        : target_function_(element_type()) {};

      /**
       * \brief Construct wrapper from callable
       *
       * \param target_function Callable object to wrap
       */
      explicit
      func_plusfunc_omx(element_type target_function)
        : target_function_(target_function) {};

      /**
       * \brief Evaluate the wrapped callable with the first argument x passed
       *        through the logarithm and the result divided by x
       *
       * \param x Argument to be passed through the logarithm before being
       *        passed to the wrapped callable
       * \param rest Function parameter pack that is passed on to the wrapped
       *        callable unchanged
       *
       * \return The value returned by the wrapped callable, divided by x
       */
      numeric_type operator()(numeric_type x, Trest... rest) const
      {
        using std::log;
        return(target_function_(log(static_cast<numeric_type>(1)-x), rest...)
               / (static_cast<numeric_type>(1)-x));
      };

      /**
       * \brief Evaluate the integral over the plus function
       *
       * \details
       * For Mellin convolutions, we need an extra term which corresponds to
       * \f[
       *   I(x) = \int_0^x \mathrm{d}y g(y,\dots)
       * \f]
       * Since this wrapper is supposed to model \f$1-x\f$ plus function kernels
       * and the coeffients \f$c_i(\dots)\f$ do not depend on \f$x\f$, we can
       * analytically calculate the integral and just evaluate the result:
       * \f[
       *   I(x) = \sum_k c_i(\dots) \int_0^x \mathrm{d}y \frac{\log^k(1-y)}{1-y}
       *        = \sum_k c_i(\dots) \frac{-\log^{k+1}(1-x)}{k+1}
       * \f]
       * The implementation is only valid if the wrapped polynomial has no
       * negative powers.
       *
       * \param x Upper integration bound \f$x\f$
       * \param rest Function parameter pack that is passed on to the wrapped
       *        callable unchanged
       *
       * \return The value of the integral
       */
      numeric_type eval_plus_int(numeric_type x, Trest...) const
      {
        using std::log;
        using std::pow;
        int min_power = target_function_.min_power();
        assert((min_power >= 0, "Only non-negative exponents are supported"));

        numeric_type log_omx = log(static_cast<numeric_type>(1)-x);

        // Calculate the integral over the plus function
        size_t num_coeffs = (target_function_.max_power()+1) - min_power;
        std::vector<numeric_type> plusfunc_ints(num_coeffs,
                                                static_cast<numeric_type>(0));
        // Compute the exponents that appear in the integral over the plus
        // functions
        std::iota(
          plusfunc_ints.begin(),
          plusfunc_ints.end(),
          static_cast<numeric_type>(min_power+1)
        );
        // Evaluate the integrals over the plus functions
        std::transform(
          plusfunc_ints.begin(),
          plusfunc_ints.end(),
          plusfunc_ints.begin(),
        [log_omx](numeric_type e) { return(-pow(log_omx,e)/e); }
        );

        // Combine the integrals over the individual plus functions with the
        // coefficients
        return(target_function_.eval_subst(plusfunc_ints));
      };

    private:
      element_type target_function_;
    };
  }
}

#endif
