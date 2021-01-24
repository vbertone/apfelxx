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
   * @brief The OperatorGPD class inherits from the Operator class,
   * used for "standard" Mellin convolutions, and implements the
   * operation required to compute GPD-like convolutions.
   */
  class OperatorGPD: public Operator
  {
  public:
    /**
     * @brief The OperatorGPD constructor.
     * @param gr: the Grid object
     * @param expr: the expression to be transformed
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     */
    OperatorGPD(Grid const& gr, Expression const& expr, double const& eps = 1e-5);

    /**
     * @brief Function that computes the actual operator on the grid.
     */
    void ComputeOperator() override;

    /**
     * @name Binary operators
     */
    ///@{
    Distribution operator *= (Distribution const& d) const override; //!< this *= Distribution
    OperatorGPD& operator *= (OperatorGPD const& o);                 //!< this *= Operator
    ///@}
  };
}
