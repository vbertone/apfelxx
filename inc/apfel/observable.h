//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/set.h"
#include "apfel/operator.h"

namespace apfel
{
  /**
   * @brief The Observable class encapsulates sets of operators and
   * sets of operators for an easy computation of observebles deriving
   * from the covolution of the two.
   */
  class Observable
  {
  public:

    Observable() =  delete;

    /**
     * @brief Observable default constructor.
     * @param CoefficientFunctions: a Set<Operator>-valued function returning the operators
     * @param Distributions: a Set<Distribution>-valued function returning the distributions
     */
    Observable(std::function<Set<Operator>(double const&)>     const& CoefficientFunctions,
	       std::function<Set<Distribution>(double const&)> const& Distributions);

    /**
     * @name Functions that evaluate the the observable at the scale
     * Q.
     */
    ///@{
    /**
     * @brief This function returns the observable as a distribution.
     * @param Q: the scale where the observable has to be evaluated
     * @return the observable in Q as a distribution
     */
    Distribution Evaluate(double const& Q)                  const;
    /**
     * @brief This function returns the observable in Q interpolated in x.
     * @param x: the value to be interpolate on the x-space grid
     * @param Q: the scale where the observable has to be evaluated
     * @return the observable in Q interpolated in x
     */
    double       Evaluate(double const& x, double const& Q) const;
    ///@}

    /**
     * @brief Set the set of ditributions keeping the same set of
     * coefficient functions.
     * @param Distributions: the new set of distributions
     */
    void SetDistributions(std::function<Set<Distribution>(double const&)> const& Distributions) { _Distributions = Distributions; }

  private:
    std::function<Set<Operator>(double const&)>     _CoefficientFunctions;
    std::function<Set<Distribution>(double const&)> _Distributions;
  };

}
