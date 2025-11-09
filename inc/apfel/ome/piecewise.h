/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Callable for modelling piecewise functions
 */

#ifndef LIBOME_PIECEWISE_H
#define LIBOME_PIECEWISE_H

#include <stdexcept>
#include <vector>
#include <algorithm>

namespace apfel
{
  namespace ome
  {
    /**
     * \brief Exception: definition of piecewise function is inconsistent
     */
    class inconsistent_piecewise : public std::logic_error
    {
    public:
      // Inherit constructors from base class
      using std::logic_error::logic_error;
    };

    /**
     * \brief One-dimensional piecewise function
     *
     * \details
     * This class models a piecewise function defined on a subdivision of the
     * real axis into non-overlapping, contiguous itervals. Upon evaluation, the
     * first argument is used to look up on which interval the function is
     * evaluated and the corresponding function, stored as a callable is
     * evaluated.
     *
     * \tparam Tnum Numerical type for the first argument upon evaluation, the
     *         return type
     * \tparam Tfunc Callable type for the functions on each interval
     * \tparam Trest Pack of types for the remaining arguments of the wrapped
     *         callable
     */
    template<typename Tnum, typename Tfunc, typename... Trest>
    class piecewise
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
       * Constructs a function with only one interval and default constructed
       * callable.
       */
      piecewise()
        : boundaries_(), pieces_{{element_type()}} {};

      /**
       * \brief Construct piecewise function from interval boundaries and
       *        callables on those intervals.
       *
       * \details
       * The vector boundaries must be sorted in ascending order and have
       * exactly one element less that the vector pieces. The first interval
       * starts at negative infinity and extends up to, but excluding the first
       * boundary \f$b_0\f$. Each further interval is given by
       * \f$[b_i,b_{i+1})\f$ and the last interval extends to positive infinity.
       *
       * \param boundaries Vector of interval boundaries
       * \param pieces Vector of callable objects representing the function
       *        on each interval
       */
      piecewise(const std::vector<numeric_type> boundaries, const std::vector<element_type> pieces)
        : boundaries_(boundaries), pieces_(pieces)
      {
        if(boundaries_.size()+1 != pieces_.size())
          throw(inconsistent_piecewise("The must be one more function than boundaries in a piecewise function"));

        if(!std::is_sorted(boundaries_.begin(), boundaries_.end()))
          throw(inconsistent_piecewise("Interval boundaries must be specified in ascending order"));
      };

      /**
       * \brief Helper to allow construction directly from intializer lists
       *
       * \details
       * See \ref piecewise(const std::vector<numeric_type>,
       * const std::vector<element_type>) for details.
       */
      piecewise(std::initializer_list<numeric_type> boundaries, std::initializer_list<element_type> pieces)
        : piecewise(std::vector<numeric_type>(boundaries), std::vector<element_type>(pieces))
      {};

      /**
       * \brief Evaluate piecewise function
       *
       * \param x First argument that decides on which interval we have to evaluate
       * \param rest Function parameter pack that is passed on to the wrapped
       *        callable unchanged
       *
       * \return The value of the piecewise function at the given point
       */
      numeric_type operator()(numeric_type x, Trest... rest) const
      {
        if(pieces_.empty())
          return(static_cast<numeric_type>(0));

        auto piece = pieces_.cbegin();
        for(auto boundary = boundaries_.cbegin(); boundary != boundaries_.cend(); ++piece, ++boundary)
          {
            if(x < *boundary)
              {
                return((*piece)(x, rest...));
              }
          }
        return((*piece)(x, rest...));
      };

    private:
      std::vector<numeric_type> boundaries_;
      std::vector<element_type> pieces_;
    };
  };
}
#endif
