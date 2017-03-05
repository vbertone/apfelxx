//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"

using namespace std;

namespace apfel
{
  /**
   * @brief Identity expression (delta function)
   */
  class Identity: public Expression
  {
  public:
  Identity(): Expression() { }
    double Local(double const& x) const { return 1 + 0 * x; }
  };

  /**
   * @brief Zero expression
   */
  class Null: public Expression
  {
  public:
  Null(): Expression() { }
  };

}
