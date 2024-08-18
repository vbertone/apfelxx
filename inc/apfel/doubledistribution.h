//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/distribution.h"
#include "apfel/matrix.h"
#include "apfel/lagrangeinterpolator.h"

namespace apfel
{
  /**
   * @brief The DoubleDistribution class defines one of the basic
   * objects of APFEL++. This is essentially the discretisation of a
   * two-dimensional function that can be conveniently used for
   * convolutions.
   */
  class DoubleDistribution
  {
  public:
    DoubleDistribution() = delete;

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DoubleDistribution constructor
     * @param g1: the Grid object that defines the first x-space interpolation grid
     * @param g2: the Grid object that defines the second x-space interpolation grid
     * @param InDistFunc: a function of x1 and x2 to be tabulated
     */
    DoubleDistribution(Grid                                                const& g1,
                       Grid                                                const& g2,
                       std::function<double(double const&, double const&)> const& InDistFunc);

    /**
     * @brief The DoubleDistribution copy constructor
     * @param obj: object to be copied
     */
    DoubleDistribution(DoubleDistribution const& obj);

    /**
     * @brief The DoubleDistribution constructor
     * @param d1: first single-distribution
     * @param d2: second single-distribution
     */
    DoubleDistribution(Distribution const& d1, Distribution const& d2);
    ///@}

    /**
     * @name Evaluate functions
     * List of functions that perform the interpolation on the x-space
     * grids.
     */
    ///@{
    /**
     * @brief Function that evaluates the interpolated function on the joint grid.
     * @param x1: the value in x1 to be interpolated
     * @param x2: the value in x2 to be interpolated
     * @return the interpolated result
     */
    double Evaluate(double const& x1, double const& x2) const;

    /**
     * @brief Function that evaluates the derivative of the
     * interpolated function on the joint grid.
     * @param x1: the value in x1 where the derivative has to be computed
     * @param x2: the value in x2 where the derivative has to be computed
     * @return the derivative of the interpolated function
     */
    double Derive(double const& x1, double const& x2) const;

    /**
     * @brief Function that evaluates the integral of the interpolated
     * function in the interval x1 in [a1,b1] and x2 in [a2,b2] on the
     * joint grid.
     * @param a1: the lower integration bound for the first variable
     * @param b1: the upper integration bound for the first variable
     * @param a2: the lower integration bound for the second variable
     * @param b2: the upper integration bound for the second variable
     * @return the integral of the interpolated function
     */
    double Integrate(double const& a1, double const& b1, double const& a2, double const& b2) const;
    ///@}

    /**
     * @name Getters
     */
    ///@{
    Grid                                     const& GetFirstGrid()  const { return _g1; }
    Grid                                     const& GetSecondGrid() const { return _g2; }
    std::vector<std::vector<matrix<double>>> const  GetDistributionSubGrid()   const { return _dDSubGrid; }
    matrix<double>                           const  GetDistributionJointGrid() const { return _dDJointGrid; }
    ///@}

  private:
    Grid                              const& _g1;           //!< The first interpolation grid
    Grid                              const& _g2;           //!< The second interpolation grid
    LagrangeInterpolator              const  _li1;          //!< The first Lagrange interpolator
    LagrangeInterpolator              const  _li2;          //!< The second Lagrange interpolator
    std::vector<std::vector<matrix<double>>> _dDSubGrid;    //!< The array with the double distribution values on the subgrids
    matrix<double>                           _dDJointGrid;  //!< The array with the double distribution values on the joint grids

    friend std::ostream& operator << (std::ostream& os, DoubleDistribution const& sg);
  };

  /**
   * @brief Method which prints Interpolator with cout <<.
   */
  std::ostream& operator << (std::ostream& os, DoubleDistribution const& in);
}
