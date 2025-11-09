/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Evaluatable (and nestable) Laurent polynomials
 */

#ifndef LIBOME_LAURENT_POLYNOMIAL_H
#define LIBOME_LAURENT_POLYNOMIAL_H

#include <utility>
#include <type_traits>
#include <initializer_list>
#include <iterator>
#include <vector>
#include <cmath>
#include <numeric>
#include <apfel/ome/traits.h>

namespace apfel
{
  namespace ome
  {
    /**
     * \brief View on an evaluatable Laurent polynomial
     *
     * \details
     * This class provides a view on a Laurent polyonmial (i.e., it does not store
     * the Laurent polynomial itself, but it only refers to one). It can be nested
     * (see documentation of \ref laurent_polynomial) and it can be evaluated.
     *
     * \tparam Tnum Numerical data type. This is the type of the polynomial
     *         variable \f$x\f$, the type of the coefficients or the return type
     *         of evaluating the coefficients \f$c_i(\dots)\f$ and the type that
     *         the polynomial evaluates to.
     * \tparam Tcoeff Type of the coefficients \f$c_i\f$ of the polynomial. This
     *         can either be a numerical type that can be implicitly converted to
     *         Tnum or a callable that evaluates to Tnum when invoked.
     * \tparam Tcoefffargs Pack of types that specifies the types of the arguments
     *         when invoking Tcoeff. If Tcoeff is a numerical type (not a
     *         callable) this pack should be empty.
     */
    template<typename Tnum, typename Tcoeff, typename... Tcoeffargs>
    class laurent_polynomial_view
    {
    public:
      /// Type alias for the numerical type template parameter
      using numeric_type = Tnum;
      /// Type alias for the coefficient type template parameter
      using coefficient_type = Tcoeff;
      /// Boolean type alias for whether the coefficient_type has a view
      using coefficient_has_view = has_view<coefficient_type>;
      /// Boolean type alias indicating that this class has an eval_plus_int method
      using has_eval_plus_int = std::true_type;
      /// Type alias for the view of the coefficient type (if it exists, otherwise it is void)
      using coefficient_view_type = typename get_view_type<coefficient_type>::type;
      /// Type alias for the coefficient iterators
      using const_reverse_iterator = typename std::vector<coefficient_type>::const_reverse_iterator;

      /**
       * \brief Construct from const iterators
       *
       * \param crbegin Reverse iterator to last coefficient
       * \param crend Reverse iterator to behind first coefficient
       * \param min_power Exponent of the lowest-order monomial (\f$n_0\f$).
       */
      laurent_polynomial_view(
        const_reverse_iterator crbegin,
        const_reverse_iterator crend,
        int min_power = 0)
        : coeffs_crbegin_(crbegin),
          coeffs_crend_(crend),
          min_power_(min_power) {};

      /**
       * \brief Swap function
       */
      friend
      void swap(laurent_polynomial_view& first, laurent_polynomial_view& second)
      {
        using std::swap;
        swap(first.coeffs_crbegin_,second.coeffs_crbegin_);
        swap(first.coeffs_crend_,second.coeffs_crend_);
        swap(first.min_power_,second.min_power_);
      };

      /**
       * \brief Evaluate Laurent polynomial
       *
       * \details
       * Evaluate the polynomial for a given set of values of the variables. The
       * first argument is the value for the variable \f$x\f$. If there are more
       * arguments, they are passed on to the evaluation of each coefficient.
       *
       * \param x Value for the polynomial variable \f$x\f$.
       * \param rest Pack of function arguments that are passed to the
       *        evaluation of the coefficients \f$c_i(\dots)\f$.
       *
       * \return Numerical value for the Laurent polynomial evaluated at the
       *         given point.
       */
      numeric_type operator()(numeric_type x, Tcoeffargs... rest) const
      {
        using std::pow;
        const numeric_type prefactor = (min_power_ != 0) ? pow(x, min_power_) : static_cast<numeric_type>(1);

        if constexpr(sizeof...(rest) > 0)
          {
            if(std::next(coeffs_crbegin_) == coeffs_crend_)
              {
                return(prefactor*(*coeffs_crbegin_)(rest...));
              }
            else
              {
                return(prefactor*std::accumulate(
                         coeffs_crbegin_,
                         coeffs_crend_,
                         static_cast<numeric_type>(0),
                         [&x, &rest...] (const numeric_type& acc, const coefficient_type& next)
                {
                  return(x*acc + next(rest...));
                }
                       ));
              }
          }
        else
          {
            if(std::next(coeffs_crbegin_) == coeffs_crend_)
              {
                return(prefactor*(*coeffs_crbegin_));
              }
            else
              {
                return(prefactor*std::accumulate(
                         coeffs_crbegin_,
                         coeffs_crend_,
                         static_cast<numeric_type>(0),
                         [&x] (const numeric_type& acc, const coefficient_type& next)
                {
                  return(x*acc + next);
                }
                       ));
              }
          }
      };

      /**
       * \brief Evaluate polynomial with precomputed monomials
       *
       * \details
       * Evaluate the Laurent polynomial, but instead of the \f$x^i\f$, the
       * monomials are taken from the input vector subst. Effectively this method computes the inner product between subst and the vector of polynomial coefficients.
       *
       * \param subst Vector of precomputed monomials. The caller has to
       *        ensure that it has the same length as the vector of
       *        polynomial coefficients
       * \param rest Pack of function arguments that are passed to the
       *        evaluation of the coefficients \f$c_i(\dots)\f$.
       *
       * \return Numerical value for the Laurent polynomial evaluated with the
       *         given monomials.
       */
      numeric_type eval_subst(const std::vector<numeric_type>& subst, Tcoeffargs... rest) const
      {
        if constexpr(sizeof...(rest) > 0)
          {
            return(std::transform_reduce(
                     coeffs_crbegin_,
                     coeffs_crend_,
                     subst.crbegin(),
                     static_cast<numeric_type>(0),
                     std::plus<>(),
                     [&rest...] (const coefficient_type& c, const numeric_type& s)
            {
              return(c(rest...) * s);
            }
                   ));
          }
        else
          {
            return(std::transform_reduce(
                     coeffs_crbegin_,
                     coeffs_crend_,
                     subst.crbegin(),
                     static_cast<numeric_type>(0)
                   ));
          }
      };

      /**
       * \brief Evaluate polynomial, but call eval_plus_int on wrapped type
       *
       * \details
       * This is meant to be used in conjunction with func_plusfunc_omx, which
       * allows to calculate the "extra" integral of plus functions over
       * \f$[0,x]\f$ analytically. This method exposes that functionality if
       * func_plusfunc_omx is wrapped into further laurent_polynomials.
       *
       * \param x Value for the polynomial variable \f$x\f$.
       * \param rest Pack of function arguments that are passed to the
       *        evaluation of the coefficients \f$c_i(\dots)\f$.
       *
       * \return Numerical value for the Laurent polynomial evaluated for the
       *         given values of the variables, but with eval_plus_int called
       *         on the wrapped callable.
       */
      numeric_type eval_plus_int(numeric_type x, Tcoeffargs... rest) const
      {
        using std::pow;
        const numeric_type prefactor = (min_power_ != 0) ? pow(x, min_power_) : static_cast<numeric_type>(1);

        if(std::next(coeffs_crbegin_) == coeffs_crend_)
          {
            return(prefactor * coeffs_crbegin_->eval_plus_int(rest...));
          }
        else
          {
            return(prefactor*std::accumulate(
                     coeffs_crbegin_,
                     coeffs_crend_,
                     static_cast<numeric_type>(0),
                     [&x, &rest...] (const numeric_type& acc, const coefficient_type& next)
            {
              return(x*acc + next.eval_plus_int(rest...));
            }
                   ));
          }
      };

      /**
       * \brief Truncate the Laurent polynomial to a given order
       *
       * \details
       * Truncation means that the polynomial is only evaluated up to a certain
       * order, which can by lower than the highest available term. If the
       * requested truncation order is lower than the minimal exponent
       * (\f$n_0\f$), the resulting truncated polynomial will always evaluate to
       * zero. If the requested truncation order is higher than the maximal
       * exponent, no truncation will be performed.
       *
       * \param truncation_order Exponent of the last order to include
       *
       * \return A laurent_polynomial_view of the truncated polynomial
       */
      laurent_polynomial_view truncate(int truncation_order) const
      {
        const int requested_orders = truncation_order - min_power_ + 1;
        const int available_orders = coeffs_crend_ - coeffs_crbegin_;

        if(requested_orders <= 0)
          {
            return(laurent_polynomial_view(coeffs_crbegin_,
                                           coeffs_crbegin_, min_power_));
          }
        else if(requested_orders > available_orders)
          {
            return(laurent_polynomial_view(coeffs_crbegin_,
                                           coeffs_crend_, min_power_));
          }
        else
          {
            return(laurent_polynomial_view(
                     coeffs_crbegin_ + (available_orders - requested_orders),
                     coeffs_crend_,
                     min_power_
                   ));
          }
      };

      /**
       * \brief Clone the view on this Laurent polynomial
       *
       * \return A laurent_polynomial_view on this polynomial
       */
      laurent_polynomial_view get_view() const
      {
        return(laurent_polynomial_view(*this));
      };

      /**
       * \brief Get the minimum exponent \f$n_0\f$ of the polynomial
       *
       * \return Minimum exponent
       */
      int min_power() const
      {
        return(min_power_);
      };

      /**
       * \brief Get the maximum exponent \f$N\f$ of the polynomial
       *
       * \details
       * If the polynomial is empty (i.e. it doesn't have any coefficients)
       * it returns min_power - 1
       *
       * \return Maximum exponent
       */
      int max_power() const
      {
        return(min_power_ + (coeffs_crend_ - coeffs_crbegin_ - 1));
      };

      /**
       * \brief Access to the coefficients of the Laurent polynomial
       *
       * \param pow Index of the coefficient to access
       *
       * \return Constant reference to the coefficient of \f$x^\mathtt{pow}\f$.
       */
      const coefficient_type& operator[](int pow) const
      {
        const int requested_position = pow - min_power_;
        const int size = coeffs_crend_ - coeffs_crbegin_;
        if(requested_position < 0 || requested_position >= size)
          return zero_value_;
        else
          return *(coeffs_crbegin_ + ((size-1) - requested_position));
      };

    protected:
      /// Reverse iterator pointing to the beginning of the coefficient vector
      const_reverse_iterator coeffs_crbegin_;
      /// Reverse iterator pointing to the end of the coefficient vector
      const_reverse_iterator coeffs_crend_;

      /**
       * \brief Default constructor
       */
      laurent_polynomial_view(int min_power = 0)
        : min_power_(min_power) {};

    private:
      int min_power_;

      inline const static coefficient_type zero_value_{};
    };

    /**
     * \brief Potentially nested, evaluatable Laurent polynomial
     *
     * \details
     * This class implements a Laurent polynomial \f$p(x, \dots)\f$ (i.e. it can
     * have monomials with negative degree).
     * \f[
     *    p(x, \dots) = \sum_{n=n_0}^N x^i c_i(\dots)
     * \f]
     * The polynomial can be nested, which means that its coefficients \f$c_i\f$
     * can be callables. In particular, they can be laurent_polynomials. This
     * allows to implement multivariate polynomials. Coefficients are specified
     * at the time of construction and cannot be modified.
     *
     * \tparam Tnum Numerical data type. This is the type of the polynomial
     *         variable \f$x\f$, the type of the coefficients or the return type
     *         of evaluating the coefficients \f$c_i(\dots)\f$ and the type that
     *         the polynomial evaluates to.
     * \tparam Tcoeff Type of the coefficients \f$c_i\f$ of the polynomial. This
     *         can either be a numerical type that can be implicitly converted to
     *         Tnum or a callable that evaluates to Tnum when invoked.
     * \tparam Tcoefffargs Pack of types that specifies the types of the arguments
     *         when invoking Tcoeff. If Tcoeff is a numerical type (not a
     *         callable) this pack should be empty.
     */
    template<typename Tnum, typename Tcoeff, typename... Tcoeffargs>
    class laurent_polynomial
      : public laurent_polynomial_view<Tnum, Tcoeff, Tcoeffargs...>
    {
    public:
      /// Type alias for the underlying view from which we inherit
      using view_type = laurent_polynomial_view<Tnum, Tcoeff, Tcoeffargs...>;
      // Make type aliases from base class usable here
      using typename view_type::numeric_type;
      using typename view_type::coefficient_type;


      /**
       * \brief Construct a null polynomial (always evaluates to zero)
       */
      laurent_polynomial()
        : view_type(0),
          coefficients_()
      {
        this->coeffs_crbegin_ = coefficients_.crbegin();
        this->coeffs_crend_ = coefficients_.crend();
      };

      /**
       * \brief Construct Laurent polynomial from a vector of coefficients
       *
       * \param coefficients Vector of coefficients \f$c_i(\dots)\f$ of the
       *        polynomial. The first element is the coefficient of the
       *        lowest-order monomial \f$x^{n_0}\f$ and subsequent elements are
       *        the coefficients of the higher powers.
       * \param min_power Exponent of the lowest-order monomial (\f$n_0\f$).
       */
      laurent_polynomial(const std::vector<coefficient_type>& coefficients,
                         int min_power = 0)
        : view_type(min_power),
          coefficients_(coefficients)
      {
        this->coeffs_crbegin_ = coefficients_.crbegin();
        this->coeffs_crend_ = coefficients_.crend();
      };

      /**
       * \brief Helper to construct a Laurent polynomial directly from an
       *        initalizer list.
       *
       * \details
       * See \ref laurent_polynomial(const std::vector<coefficient_type>&, int) for details.
       */
      laurent_polynomial(std::initializer_list<coefficient_type> coefficients,
                         int min_power = 0)
        : view_type(min_power),
          coefficients_(coefficients)
      {
        this->coeffs_crbegin_ = coefficients_.crbegin();
        this->coeffs_crend_ = coefficients_.crend();
      };

      /**
       * \brief Copy constructor
       */
      laurent_polynomial(const laurent_polynomial& other)
        : view_type(other),
          coefficients_(other.coefficients_)
      {
        this->coeffs_crbegin_ = coefficients_.crbegin();
        this->coeffs_crend_ = coefficients_.crend();
      };

      /**
       * \brief Swap function
       */
      friend
      void swap(laurent_polynomial& first, laurent_polynomial& second)
      {
        // enable proper ADL
        using std::swap;

        // swap base class members
        swap(static_cast<view_type&>(first),static_cast<view_type&>(second));

        // swap local data members
        swap(first.coefficients_,second.coefficients_);

        // update iterators
        first.coeffs_crbegin_ = first.coefficients_.crbegin();
        first.coeffs_crend_   = first.coefficients_.crend();
        second.coeffs_crbegin_ = second.coefficients_.crbegin();
        second.coeffs_crend_   = second.coefficients_.crend();
      };

      /**
       * \brief Assignment operator
       */
      laurent_polynomial& operator=(laurent_polynomial other)
      {
        swap(*this, other);
        return(*this);
      };

      /**
       * \brief Move operator
       */
      laurent_polynomial(laurent_polynomial&& other)
        : laurent_polynomial()
      {
        swap(*this, other);
      };

    private:
      std::vector<coefficient_type> coefficients_;
    };

  }
}

#endif
