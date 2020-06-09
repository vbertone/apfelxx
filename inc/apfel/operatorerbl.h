//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/operator.h"

namespace apfel
{
  /**
   * @brief The class, derived from the Operator class, defines a
   * special kind of Object necessary for the evolution of GPDs in the
   * ERBL region.
   */
  class OperatorERBL: public Operator
  {
  public:

    OperatorERBL() = delete;

    /**
     * @brief The Operator constructor.
     * @param gr: the Grid object
     * @param expr: the expression to be transformed
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     */
    OperatorERBL(Grid const& gr, Expression const& expr, double const& eps = 1e-5);
  };
}
