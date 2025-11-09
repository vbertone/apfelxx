/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Container for mixed distributions with regular, plus and delta part
 */

#ifndef LIBOME_RPD_DISTRIBUTION_H
#define LIBOME_RPD_DISTRIBUTION_H

#include <optional>

namespace apfel
{
  namespace ome
  {
    /**
     * \brief Mixed distribution with regular, plus and delta part
     *
     * \details
     * This is a container that models a distribution \f$f(\dots,x)\f$ which is a
     * sum of a regular, a plus and a delta part:
     * \f[
     *   f(\dots,x) = f_\mathrm{reg}(\dots,x)
     *     +\left[f_+(\dots,x)\right]_+
     *     +\delta(1-x) f_\delta(\dots)
     * \f]
     * Whether or not each of the three parts is present is optional. All three
     * parts should have the compatible evaluation signatures: the parameters
     * denoted by \f$\dots\f$ above should be the same for each of the three
     * parts.
     *
     * \tparam Tfuncreg Type for the wrapped regular part
     * \tparam Tfuncplus Type for the wrapped plus part
     * \tparam Tfuncdelta Type for the wrapped delta part
     */
    template<typename Tfuncreg,
             typename Tfuncplus,
             typename Tfuncdelta>
    class rpd_distribution
    {
    public:
      /// Type alias for the regular part template parameter
      using element_regular_type = Tfuncreg;
      /// Type alias for the plus part template parameter
      using element_plus_type = Tfuncplus;
      /// Type alias for the delta part template parameter
      using element_delta_type = Tfuncdelta;


      /**
       * \brief Default construtor
       *
       * \details
       * Constructs an "empty" container view with no regular, plus or delta
       * part.
       */
      rpd_distribution()
        : regular_(std::nullopt), plus_(std::nullopt), delta_(std::nullopt) {};

      /**
       * \brief Construct distribution from regular, plus and delta parts
       *
       * \details
       * The individual parts are wrapped in std::optional so that their absence
       * can be signalled by passing std::nullopt.
       *
       * \param regular Regular part. Can be a numerical value or a callable
       * \param plus Plus part. Can be a numerical value or a callable
       * \param delta Delta part. Can be a numerical value or a callable
       */
      rpd_distribution(
        std::optional<element_regular_type> regular,
        std::optional<element_plus_type> plus,
        std::optional<element_delta_type> delta
      )
        : regular_(regular), plus_(plus), delta_(delta) {};

      /**
       * \brief Check whether the regular part is present
       *
       * \return True if the regular part is present, false otherwise.
       */
      bool has_regular() const
      {
        return(regular_.has_value());
      }

      /**
       * \brief Check whether the plus part is present
       *
       * \return True if the plus part is present, false otherwise.
       */
      bool has_plus() const
      {
        return(plus_.has_value());
      }

      /**
       * \brief Check whether the delta part is present
       *
       * \return True if the delta part is present, false otherwise.
       */
      bool has_delta() const
      {
        return(delta_.has_value());
      }

      /**
       * \brief Retrieve the regular part
       *
       * \return Constant reference to an std::optional wrapping the regular
       *         part. If the regular part is absent, the std::optional
       *         contains std::nullopt.
       */
      const std::optional<element_regular_type>& get_regular() const
      {
        return(regular_);
      };

      /**
       * \brief Retrieve the plus part
       *
       * \return Constant reference to an std::optional wrapping the plus part.
       *         If the plus part is absent, the std::optional contains
       *         std::nullopt.
       */
      const std::optional<element_plus_type>& get_plus() const
      {
        return(plus_);
      };

      /**
       * \brief Retrieve the delta part
       *
       * \return Constant reference to an std::optional wrapping the delta part.
       *         If the delta part is absent, the std::optional contains
       *         std::nullopt.
       */
      const std::optional<element_delta_type>& get_delta() const
      {
        return(delta_);
      };

      /**
       * \brief Extract the n-th coefficient from the regular, plus and delta parts
       *        into a new rpd_distribution container
       *
       * \details
       * In order to be applicable, the regular, plus and delta parts must have a
       * coefficient_type, coefficient_has_view and coefficient_view_type member
       * types. Also the n-th coefficient must be extractable using operator[][n].
       * If coefficient_has_view is std::true_type, the function extracts the view
       * of the n-th coefficient instead of the coefficient itself.
       *
       * \param n Exponent of the coefficient to extract
       *
       * \return rpd_distribution container with the n-th coefficient (or a view
       *         on it, if available) of the regular, plus and delta parts
       */
      auto get_coefficient(int n) const
      {
        using reg_type
        = std::conditional_t<element_regular_type::coefficient_has_view::value,
        typename element_regular_type::coefficient_view_type,
        typename element_regular_type::coefficient_type>;
        using plus_type
        = std::conditional_t<element_plus_type::coefficient_has_view::value,
        typename element_plus_type::coefficient_view_type,
        typename element_plus_type::coefficient_type>;
        using delta_type
        = std::conditional_t<element_delta_type::coefficient_has_view::value,
        typename element_delta_type::coefficient_view_type,
        typename element_delta_type::coefficient_type>;
        std::optional<reg_type> reg;
        std::optional<plus_type> plus;
        std::optional<delta_type> delta;

        if(has_regular())
          {
            if constexpr (element_regular_type::coefficient_has_view::value)
              reg = std::make_optional((*regular_)[n].get_view());
            else
              reg = std::make_optional((*regular_)[n]);
          }
        if(has_plus())
          {
            if constexpr (element_plus_type::coefficient_has_view::value)
              plus = std::make_optional((*plus_)[n].get_view());
            else
              plus = std::make_optional((*plus_)[n]);
          }
        if(has_delta())
          {
            if constexpr (element_delta_type::coefficient_has_view::value)
              delta = std::make_optional((*delta_)[n].get_view());
            else
              delta = std::make_optional((*delta_)[n]);
          }

        return(
                rpd_distribution<reg_type, plus_type, delta_type> {reg, plus, delta}
              );
      };

      /**
       * \brief Apply the truncate operation on each of the regular, plus and
       *        delta parts (if supported)
       *
       * \param n Exponent of the last order to include
       *
       * \return An rpd_distribution with truncated regular, plus and delta
       *         parts
       */
      auto truncate(int n) const
      {
        using reg_type = decltype(std::declval<element_regular_type>().truncate(n));
        using plus_type = decltype(std::declval<element_plus_type>().truncate(n));
        using delta_type = decltype(std::declval<element_delta_type>().truncate(n));

        std::optional<reg_type> reg;
        std::optional<plus_type> plus;
        std::optional<delta_type> delta;

        if(has_regular())
          {
            reg = std::make_optional(regular_->truncate(n));
          }
        if(has_plus())
          {
            plus = std::make_optional(plus_->truncate(n));
          }
        if(has_delta())
          {
            delta = std::make_optional(delta_->truncate(n));
          }

        return(
                rpd_distribution<reg_type, plus_type, delta_type> {reg, plus, delta}
              );
      };

    private:
      std::optional<element_regular_type> regular_;
      std::optional<element_plus_type> plus_;
      std::optional<element_delta_type> delta_;
    };

  }
}

#endif
