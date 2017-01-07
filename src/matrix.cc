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
  matrix<T>::matrix():
    _size{0,0}
  {
  }

  //_________________________________________________________________________
  template<class T>
  void matrix<T>::resize(int const& row, int const& col, double const& v)
  {
    _size = {row,col};
    _data.resize(row*col, v);
  }

  // type constrain
  template class matrix<int>;
  template class matrix<float>;
  template class matrix<double>;
}
