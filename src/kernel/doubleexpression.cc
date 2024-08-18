//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubleexpression.h"

namespace apfel
{
  //_________________________________________________________________________
  DoubleExpression::DoubleExpression(Expression const& expr1, Expression const& expr2):
    _expr1(expr1),
    _expr2(expr2)
  {
  }
}
