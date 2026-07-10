//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/tmcoperator.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  TMCoperator::TMCoperator(double const& M, double const& b, double const& x):
    Expression(),
    _M(M),
    _b(b),
    _x(x)
  {
  }
  double TMCoperator::Regular(double const& y) const
  {
    return - _M * _b * _x  / 2 * sqrt(y / ( 1 - y )) * j1(_M * _b * _x * sqrt(( 1 - y ) / y));
  }
  double TMCoperator::Local(double const&) const
  {
    return 1;
  }
}
