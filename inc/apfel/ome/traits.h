/******************************************************************************
 * This file is part of libome                                                *
 * Copyright (C) 2025 Arnd Behring, Kay Schoenwald                            *
 * SPDX-License-Identifier: GPL-3.0-or-later                                  *
 ******************************************************************************/

/**
 * \file
 * \brief Type traits for template metaprogramming
 */

#ifndef LIBOME_TRAITS_H
#define LIBOME_TRAITS_H

#include <type_traits>

namespace apfel
{
  namespace ome
  {
    /**
     * \brief Boolean trait to check if a type has a get_view() method
     * \details Base case: false
     */
    template<typename T, typename = void>
    struct has_view : std::false_type {};

    /**
     * \brief Boolean trait to check if a type has a get_view() method
     * \details True case
     */
    template<typename T>
struct has_view<T, std::void_t<decltype(std::declval<T>().get_view())>> : std::true_type {};

    /**
     * \brief Type trait to extract the return type of the get_view() method
     * \details Base case: void (no view type found)
     */
    template<typename T, typename = void>
    struct get_view_type
    {
      /// Type alias for the view type (void if no view type was found)
      using type = void;
    };

    /**
     * \brief Type trait to extract the return type of the get_view() method
     * \details Extract actual view type
     */
    template<typename T>
    struct get_view_type<T, std::void_t<decltype(std::declval<T>().get_view())>>
    {
      /// Type alias for the view type (void if no view type was found)
      using type = decltype(std::declval<T>().get_view());
    };

    /**
     * \brief Boolean trait to check if a type has the type alias has_eval_plus_int and
     *        whether it is "true" (of type std::true_type)
     * \details Base case: false
     */
    template<typename T, typename = void>
    struct has_eval_plus_int : std::false_type {};

    /**
     * \brief Boolean trait to check if a type has the type alias has_eval_plus_int and
     *        whether it is "true" (of type std::true_type)
     * \details True case
     */
    template<typename T>
    struct has_eval_plus_int<T, std::void_t<typename T::has_eval_plus_int>>
      : std::is_same<typename T::has_eval_plus_int, std::true_type> {};

  }
}

#endif
