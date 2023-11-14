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
   * sets of T-type objects for an easy computation of observables
   * deriving from the convolution of the two. This class can contain
   * an arbitrary number of such pairs that are separatately
   * convoluted and joint when the obeservable is computed by means of
   * the "Evaluate" function.
   */
  template<class T = Distribution>
  class Observable
  {
  public:
    /**
     * @brief This structure contains a pair of sets of coefficient
     * functions and of objects.
     */
    struct ConvolutionPair
    {
      ConvolutionPair(std::function<Set<Operator>(double const&)> const& C, std::function<Set<T>(double const&)> const& O): CoefficientFunctions(C), Objects(O) {}
      std::function<Set<Operator>(double const&)> CoefficientFunctions;
      std::function<Set<T>(double const&)>        Objects;
    };

    /**
     * @brief The Observable empty constructor.
     */
    Observable();

    /**
     * @brief The Observable constructor.
     * @param ConvPair: a vector of ConvolutionPair structures containing pairs of Set<Operator>-valued and Set<T>-valued functions
     */
    Observable(std::vector<ConvolutionPair> ConvPair);

    /**
     * @brief The Observable constructor.
     * @param CoefficientFunctions: a Set<Operator>-valued function returning the operators
     * @param Objects: a Set<T>-valued function returning the relevant object
     */
    Observable(std::function<Set<Operator>(double const&)> const& CoefficientFunctions,
               std::function<Set<T>(double const&)>        const& Objects);

    /**
     * @brief Function to add a convolution pair
     * @param CoefficientFunctions: a Set<Operator>-valued function returning the operators
     * @param Objects: a Set<T>-valued function returning the relevant object
     */
    void AddConvolutionPair(std::function<Set<Operator>(double const&)> const& CoefficientFunctions,
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
     * @param ip: index of the convolution-pair vector (default: 0)
     */
    void SetObjects(std::function<Set<T>(double const&)> const& Objects, int const& ip = 0);

    /**
     * @brief Get the set of coefficient functions.
     * @param ip: index of the convolution-pair vector (default: 0)
     * @return the set of coefficient functions.
     */
    std::function<Set<Operator>(double const&)> GetCoefficientFunctions(int const& ip = 0) const;

  private:
    std::vector<ConvolutionPair> _ConvPair;
  };
}
