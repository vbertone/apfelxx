//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/expression.h"

namespace apfel
{
  //_________________________________________________________________________
  Expression::Expression(bool const& ext, double const& eta):
    _ext(ext),
    _extvar(0),
    _eta(eta)
  {
  }

  //_________________________________________________________________________
  Expression::Expression(double const& eta):
    Expression{false, eta}
  {
  }
}
