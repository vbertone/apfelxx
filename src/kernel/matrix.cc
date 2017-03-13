//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/matrix.h"

namespace apfel {
  //_________________________________________________________________________
  template<class T>
  matrix<T>::matrix(const std::size_t &row, const std::size_t &col):
    _size{row,col}
  {
    if (row*col != 0)
      _data.resize(row*col);
  }

  //_________________________________________________________________________
  template<class T>
  void matrix<T>::resize(const std::size_t &row, const std::size_t &col, T const& v)
  {
    _size = {row,col};
    _data.resize(row*col, v);
  }

  // type constrain
  template class matrix<size_t>;
  template class matrix<int>;
  template class matrix<float>;
  template class matrix<double>;
}
