//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/expression.h"

namespace apfel
{
  //_________________________________________________________________________
  Expression::Expression(double const& eta, bool const& is_xdependent):
    _extvar(0),
    _eta(eta),
    _is_xdependent(is_xdependent)
  {
  }
}
