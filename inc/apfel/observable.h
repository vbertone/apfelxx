//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/set.h"
#include "apfel/operator.h"

namespace apfel
{
  /**
   * @brief The Observable class encapsulates sets of operators and
   * sets of T-type objects for an easy computation of obsarvebles
   * deriving from the convolution of the two.
   */
  template<class T = Distribution>
  class Observable
  {
  public:
    Observable() = delete;

    /**
     * @brief The Observable constructor.
     * @param CoefficientFunctions: a Set<Operator>-valued function returning the operators
     * @param Objects: a Set<T>-valued function returning the relevant object
     */
    Observable(std::function<Set<Operator>(double const&)> const& CoefficientFunctions,
               std::function<Set<T>(double const&)>        const& Objects);

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
    T Evaluate(double const& Q) const;

    /**
     * @brief This function returns the observable in Q interpolated in x.
     * @param x: the value to be interpolate on the x-space grid
     * @param Q: the scale where the observable has to be evaluated
     * @return the observable in Q interpolated in x
     */
    double Evaluate(double const& x, double const& Q) const;
    ///@}

    /**
     * @brief Set the set of ditributions keeping the same set of
     * coefficient functions.
     * @param Objects: the new set of objects
     */
    void SetObjects(std::function<Set<T>(double const&)> const& Objects) { _Objects = Objects; }

    /**
     * @brief Get the set of coefficient functions.
     * @return the set of coefficient functions.
     */
    std::function<Set<Operator>(double const&)> GetCoefficientFunctions() const { return _CoefficientFunctions; }

  private:
    std::function<Set<Operator>(double const&)> _CoefficientFunctions;
    std::function<Set<T>(double const&)>        _Objects;
  };
}
