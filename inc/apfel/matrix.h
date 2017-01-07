//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>
using std::vector;
using std::size_t;
using std::pair;

namespace apfel
{  
  /**
   * @brief The matrix class
   *
   * A simple implementation of 2d arrays based on
   * a continous memory allocation. Elements are accessible
   * throught the (i,j) operator.
   */
  template<typename T>
  class matrix
  {
  public:
    /**
     * @brief matrix constructor
     * @param row number of rows
     * @param col number of columns
     */
    matrix(size_t const& row = 0, size_t const& col = 0);

    /**
     * @brief Resizes object and set default value
     * @param row number of rows
     * @param col number of columns
     * @param v the default value
     */
    void resize(size_t const& row, size_t const& col, T const& v = 0);

    /**
     * @brief Returns the (row,col) size pair.
     */
    pair<size_t,size_t> const& size() const { return _size; }

    //operators
    T&       operator()(size_t const& i, size_t const& j)       { return _data[i+_size.first*j]; }
    T const& operator()(size_t const& i, size_t const& j) const { return _data[i+_size.first*j]; }

  private:
    pair<size_t,size_t> _size; //!< the dimension pair
    vector<T> _data;        //!< the data array
  };
}
