//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <array>
#include <cstddef>

namespace apfel
{
  /**
   * @brief The matrix class is a simple implementation of 2d arrays
   * based on a continous memory allocation. Elements are accessible
   * throught the (i,j) operator.
   */
  template<typename T>
  class matrix
  {
  public:
    /**
     * @brief The matrix constructor.
     * @param row: number of rows
     * @param col: number of columns
     */
    matrix(size_t const& row = 0, size_t const& col = 0);

    /**
     * @brief Function that resizes object and set default value.
     * @param row: number of rows
     * @param col: number of columns
     * @param v: the default value (default: 0)
     */
    void resize(size_t const& row, size_t const& col, T const& v = 0);

    /**
     * @brief Function that set all entries of the matrix to the input value.
     * @param v: the default value
     */
    void set(T const& v);

    /**
     * @brief Returns the (row,col) size pair.
     * @param dim: the dimension
     * @returns the number of raws and columns
     */
    size_t const& size(size_t const& dim) const { return _size[dim]; }

    /**
     * @brief Returns the vector of data.
     * @returns the vector of data
     */
    std::vector<T> const& data() const { return _data; }

    /**
     * @name Binary operators involving matrices
     */
    ///@{
    T&       operator()(size_t const& i, size_t const& j)       { return _data[i*_size[1]+j]; }
    T const& operator()(size_t const& i, size_t const& j) const { return _data[i*_size[1]+j]; }
    ///@}
  private:
    std::array<size_t, 2> _size; //!< the dimension pair
    std::vector<T>        _data; //!< the data array
  };
}
