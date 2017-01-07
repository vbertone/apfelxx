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
  template<typename type>
  class matrix
  {
  public:
    /**
     * @brief matrix constructor
     */
    matrix();

    /**
     * @brief Resizes object and set default value
     * @param row number of rows
     * @param col number of columns
     * @param v the default value
     */
    void resize(int const& row, int const& col, double const& v = 0);

    /**
     * @brief Returns the (row,col) size pair.
     */
    pair<size_t,size_t> const& size() const { return _size; }

    //operators
    type&       operator()(size_t const& i, size_t const& j)       { return _data[i+_size.first*j]; }
    type const& operator()(size_t const& i, size_t const& j) const { return _data[i+_size.first*j]; }

  private:
    pair<size_t,size_t> _size; //!< the dimension pair
    vector<type> _data;        //!< the data array
  };
}
