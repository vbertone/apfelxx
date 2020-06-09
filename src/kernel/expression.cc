//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/expression.h"

namespace apfel
{
  //_________________________________________________________________________
  Expression::Expression(double const& eta):
    _eta(eta),
    _extvar(0)
  {
  }
}
