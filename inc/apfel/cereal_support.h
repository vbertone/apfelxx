//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matrix.h"

#include "apfel/cereal/cereal.hpp"
#include "apfel/cereal/types/vector.hpp"

namespace apfel::cereal
{
  /**
   * @brief Non-intrusive cereal serializer for apfel::matrix<T>.
   *
   * Only the dimension pair and the contiguous data array are stored,
   * which is enough to rebuild the object through the public
   * (row, col, data) constructor. The serializer recurses naturally:
   * for T = matrix<double> the data vector is itself made of matrices,
   * each of which is (de)serialized by the very same functions, while
   * the matrix<double> leaves bottom out on std::vector<double> and
   * hit cereal's bulk-binary fast path.
   */
  template<class Archive, class T>
  void save(Archive& ar, apfel::matrix<T> const& m)
  {
    ar(static_cast<std::size_t>(m.size(0)),
       static_cast<std::size_t>(m.size(1)),
       m.data());
  }

  template<class Archive, class T>
  void load(Archive& ar, apfel::matrix<T>& m)
  {
    std::size_t row, col;
    std::vector<T> data;
    ar(row, col, data);
    m = apfel::matrix<T>(row, col, data);
  }
}
