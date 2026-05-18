//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @brief Operator that implements target mass corrections to TMD
   * matching functions.
   */
  class TMCoperator: public Expression
  {
  public:
    TMCoperator(double const& M, double const& b, double const& x);
    double Regular(double const&) const;
    double Local(double const&)   const;
  private:
    double const _M;
    double const _b;
    double const _x;
  };
}
