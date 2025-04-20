//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

#include <vector>

namespace apfel
{
  /**
   * @defgroup UnpSFQED Unpolarised splitting functions for QED evolution
   * @ingroup SLSplittings
   */
  ///@{
  ///@}
  /**
   * @defgroup LOunpsf LO splitting functions
   * @ingroup UnpSFQED
   */
  ///@{
  /**
   * @brief Space-like O(&alpha;) non-singlet unpolarised splitting
   * function.
   */
  class P0qedns: public Expression
  {
  public:
    P0qedns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief Space-like O(&alpha;) quark-gamma unpolarised splitting
   * function.
   */
  class P0qedqgm: public Expression
  {
  public:
    P0qedqgm();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;) gamma-quark unpolarised splitting
   * function.
   */
  class P0qedgmq: public Expression
  {
  public:
    P0qedgmq();
    double Regular(double const& x) const;
  };

  /**
   * @brief Space-like O(&alpha;) gamma-gamma unpolarised splitting
   * function.
   */
  class P0qedgmgm: public Expression
  {
  public:
    P0qedgmgm();
    double Local(double const& x) const;
  };
  ///@}
}
