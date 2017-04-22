//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"
#include "apfel/set.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"

#include <functional>

using std::function;

namespace apfel
{
  /**
   * @brief The Observable class.
   */
  class Observable
  {
  public:

    Observable() =  delete;

    /**
     * @brief Observable default constructor (ZM)
     */
    Observable(function<Set<Operator>(double const&)>     const& CoefficientFunctions,
	       function<Set<Distribution>(double const&)> const& Distributions);

    /**
     * @brief Observable default constructor (ZM)
     */
    Observable(function<Distribution(double const&)> const& Obs);

    /**
     * @brief Function that evaluates the the observable at the scale Q
     */
    Distribution Evaluate(double const& Q) const;

    double Evaluate(double const& x, double const& Q) const { return _Observable(Q).Evaluate(x); }

  private:
    function<Distribution(double const&)> _Observable;
  };

}
