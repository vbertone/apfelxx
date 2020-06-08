//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matrix.h"

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  matrix<T>::matrix(std::size_t const& row, std::size_t const& col):
    _size{{row, col}}
  {
    if (row * col != 0)
      _data.resize(row * col);
  }

  //_________________________________________________________________________
  template<class T>
  void matrix<T>::resize(std::size_t const& row, std::size_t const& col, T const& v)
  {
    _size = {{row, col}};
    _data.resize(row * col, v);
  }

  // Specialisations
  template class matrix<size_t>;
  template class matrix<int>;
  template class matrix<float>;
  template class matrix<double>;
  template class matrix<std::vector<int>>;
  template class matrix<std::vector<double>>;
}
