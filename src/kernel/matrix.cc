//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matrix.h"
#include "apfel/messages.h"

#include <functional>
#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  matrix<T>::matrix(size_t const& row, size_t const& col, std::vector<T> const& d):
    // *INDENT-OFF*
    _size(std::array<size_t, 2>{row, col}),
    _data(d.empty() ? std::vector<T>(row * col) : d)
    // *INDENT-ON*
  {
  }

  //_________________________________________________________________________
  template<class T>
  matrix<T>::matrix(matrix<T> const& m):
    _size(m.size()),
    _data(m.data())
  {
  }

  //_________________________________________________________________________
  template<class T>
  void matrix<T>::resize(size_t const& row, size_t const& col, T const& v)
  {
    _size = {{row, col}};
    _data.resize(row * col, v);
  }

  //_________________________________________________________________________
  template<class T>
  void matrix<T>::set(T const& v)
  {
    std::fill(_data.begin(), _data.end(), v);
  }

  //_________________________________________________________________________
  template<class T>
  matrix<T>& matrix<T>::operator = (matrix<T> const& m)
  {
    _size = m.size();
    _data = m.data();

    return *this;
  }

  // Specialisations
  template class matrix<double>;
  template class matrix<std::vector<int>>;
  template class matrix<std::vector<double>>;
  template class matrix<matrix<double>>;

  //_________________________________________________________________________
  template<>
  matrix<double>& matrix<double>::operator += (matrix<double> const& m)
  {
    if (m.size() != _size)
      throw std::runtime_error(error("matrix<double>::operator+=", "Dimensions of matrices do not match"));

    std::transform(_data.begin(), _data.end(), m.data().begin(), _data.begin(), std::plus<double>());

    return *this;
  }

  //_________________________________________________________________________
  template<>
  matrix<double>& matrix<double>::operator -= (matrix<double> const& m)
  {
    if (m.size() != _size)
      throw std::runtime_error(error("matrix<double>::operator-=", "Dimensions of matrices do not match"));

    std::transform(_data.begin(), _data.end(), m.data().begin(), _data.begin(), std::minus<double>());

    return *this;
  }

  //_________________________________________________________________________
  template<>
  matrix<double>& matrix<double>::operator *= (double const& f)
  {
    std::transform(_data.begin(), _data.end(), _data.begin(), std::bind(std::multiplies<double>(), std::placeholders::_1, f));

    return *this;
  }

  //_________________________________________________________________________
  template<>
  matrix<double>& matrix<double>::operator /= (double const& f)
  {
    std::transform(_data.begin(), _data.end(), _data.begin(), std::bind(std::divides<double>(), std::placeholders::_1, f));

    return *this;
  }

  //_________________________________________________________________________
  template<>
  matrix<double>& matrix<double>::operator *= (matrix<double> const& m)
  {
    if (m.size(0) != _size[1])
      throw std::runtime_error(error("matrix<double>::operator*=", "Dimensions of matrices do not match"));

    // Allocate product matrix
    std::vector<double> p(_size[0] * m.size(1));

    // Compute product
    for (size_t i = 0; i < _size[0]; i++)
      for (size_t j = 0; j < m.size(1); j++)
        for (size_t k = 0; k < _size[1]; k++)
          p[i * m.size(1) + j] += _data[i * _size[1] + k] * m(k, j);

    // Adjust second dimension
    _size[1] = m.size(1);

    // Replace vector of data
    _data = p;

    return *this;
  }

  //_________________________________________________________________________
  template<>
  matrix<double> operator + (matrix<double> lhs, matrix<double> const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  template<>
  matrix<double> operator - (matrix<double> lhs, matrix<double> const& rhs)
  {
    return lhs -= rhs;
  }

  //_________________________________________________________________________
  template<>
  matrix<double> operator * (double const& s, matrix<double> rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  template<>
  matrix<double> operator * (matrix<double> lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  template<>
  matrix<double> operator / (matrix<double> lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  template<>
  matrix<double> operator * (matrix<double> lhs, matrix<double> const& rhs)
  {
    return lhs *= rhs;
  }
}
