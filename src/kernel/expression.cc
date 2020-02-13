//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/expression.h"

namespace apfel
{
  //_________________________________________________________________________
  Expression::Expression():
    _eta(1)
  {
  }

  //_________________________________________________________________________
  Expression::Expression(double const& eta):
    _eta(eta)
  {
  }
}
