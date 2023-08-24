//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <array>
#include <stddef.h>
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
     * @brief The matrix copy constructor.
     * @param m: matrix
     */
    matrix(matrix<T> const& m);

    /**
     * @brief Function that resizes object and set default value.
     * @param row: number of rows
     * @param col: number of columns
     * @param v: the default value (default: 0)
     */
    void resize(size_t const& row, size_t const& col, T const& v = 0);

    /**
     * @brief Function that sets all entries of the matrix to the
     * input value.
     * @param v: the default value
     */
    void set(T const& v);

    /**
     * @brief Returns the (row,col) size pair.
     * @param dim: the dimension
     * @returns the number of rows and columns
     */
    size_t const& size(size_t const& dim) const { return _size[dim]; }

    /**
     * @brief Returns the pair of sizes.
     * @returns the number of rows and columns
     */
    std::array<size_t, 2> const& size() const { return _size; }

    /**
     * @brief Returns the vector of data.
     * @returns the vector of data
     */
    std::vector<T> const& data() const { return _data; }

    /**
     * @name Binary operators involving matrices
     */
    ///@{
    T&       operator () (size_t const& i, size_t const& j)       { return _data[i * _size[1] + j]; }
    T const& operator () (size_t const& i, size_t const& j) const { return _data[i * _size[1] + j]; }
    matrix<T>& operator  = (matrix<T> const& m);
    matrix<T>& operator += (matrix<T> const& m);
    matrix<T>& operator -= (matrix<T> const& m);
    matrix<T>& operator *= (double const& f);
    matrix<T>& operator /= (double const& f);
    matrix<T>& operator *= (matrix<T> const& m);
    ///@}
  private:
    std::array<size_t, 2> _size; //!< The dimension pair
    std::vector<T>        _data; //!< The data array
  };

  /**
   * @name Ternary operators
   */
  ///@{
  template<class T>
  matrix<T> operator + (matrix<T> lhs, matrix<T> const& rhs); //!< matrix<T>+matrix<T>
  template<class T>
  matrix<T> operator - (matrix<T> lhs, matrix<T> const& rhs); //!< matrix<T>-matrix<T>
  template<class T>
  matrix<T> operator * (double const& s, matrix<T> rhs);      //!< Scalar*matrix<T>
  template<class T>
  matrix<T> operator * (matrix<T> lhs, double const& s);      //!< matrix<T>*Scalar
  template<class T>
  matrix<T> operator / (matrix<T> lhs, double const& s);      //!< matrix<T>/Scalar
  template<class T>
  matrix<T> operator * (matrix<T> lhs, matrix<T> const& rhs); //!< matrix<T>*matrix<T>
  ///@}
}
