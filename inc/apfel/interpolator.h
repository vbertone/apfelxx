//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"

#include <array>

namespace apfel
{
  /**
   * @brief The Interpolator class is a mother class for the x-space
   * interpolationand requires the implementation of a specialized
   * interpolation algorithm. The current version uses the joint grid
   * object stored allocated by the Grid class.
   */
  class Interpolator
  {
  public:
    virtual ~Interpolator() = default;

    /**
     * @brief The Interpolator constructor.
     * @param gr: the x-space grid object over which interpolation takes place
     */
    Interpolator(Grid const& gr);

    /**
     * @brief The Interpolator constructor.
     * @param gr: the x-space grid object over which interpolation takes place
     * @param distsubgrid: the vector of subgrids
     * @param distjointgrid: the joint subgrid
     * @note distjointgrid and distsubgrid are assumed to match the
     * structure of the grid gr.
     */
    Interpolator(Grid                             const& gr,
                 std::vector<std::vector<double>> const& distsubgrid,
                 std::vector<double>              const& distjointgrid);

    /**
     * @name Evaluate functions
     * List of functions that perform the interpolation on the x-space
     * grid. These also include the derivative and the integral of the
     * interpolated function.
     */
    ///@{
    /**
     * @brief Function that evaluates the interpolated function on the joint grid.
     * @param x: the value in x to be interpolated
     * @return the interpolated result
     */
    double Evaluate(double const& x) const;

    /**
     * @brief Function that evaluates the interpolated function on a given subgrid.
     * @param x: the value in x to be interpolated
     * @param ig: the subgrid index
     * @return the interpolated result
     */
    double Evaluate(double const& x, int const& ig) const;

    /**
     * @brief Function that evaluates the derivative of the
     * interpolated function on the joint grid.
     * @param x: the value in x where the derivative has to be computed
     * @return the derivative of the interpolated function
     */
    double Derive(double const& x) const;

    /**
     * @brief Function that evaluates the integral of the interpolated
     * function in the interval [a,b] on the joint grid.
     * @param a: the lower integration bound
     * @param b: the upper integration bound
     * @return the integral of the interpolated function
     */
    double Integrate(double const& a, double const& b) const;
    ///@}

    /**
     * @brief Pure virtual method for the interpolating functions
     * polynomial in log(x).
     * @param beta: the x-space grid index
     * @param lnx: the value (of the log) of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the interpolation weights
     * @note This interpolant is polynomial in log(x) and is used when
     * computing operators on the grid. The reason is it's translation
     * invariance on a logarithmically spaced grid that reduces the
     * ammount of computations.
     */
    virtual double InterpolantLog(int const& beta, double const& lnx, SubGrid const& sg) const = 0;

    /**
     * @brief Pure virtual method for the interpolating functions.
     * @param beta: the x-space grid index
     * @param x: the value of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the interpolation weights
     */
    virtual double Interpolant(int const& beta, double const& x, SubGrid const& sg) const = 0;

    /**
     * @brief Virtual method for the derivative of the interpolating
     * functions.
     * @param beta: the x-space grid index
     * @param lnx: the value (of the log) of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the derivarive of the interpolation weights
     */
    virtual double DerInterpolant(int const&, double const&, SubGrid const&) const { return 0; };

    /**
     * @brief Virtual method for the integral of the interpolating
     * functions.
     * @param beta: the x-space grid index
     * @param lna: the value (of the log) of the lower bound
     * @param lnb: the value (of the log) of the upper bound
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the integral of the interpolation weights
     */
    virtual double IntInterpolant(int const&, double const&, double const&, SubGrid const&) const { return 0; };

    /**
     * @brief This purely virtual function computes the lower and
     * upper bounds on which the the sum over interpolants is limited.
     * @param x: the value in x to be interpolated
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the lower and upper bounds of the grid index
     */
    virtual std::array<int,2> SumBounds(double const& x, SubGrid const& sg) const = 0;

    /**
     * @name Getters
     */
    ///@{
    Grid                             const& GetGrid()                  const { return _grid; }                  //!< The grid
    std::vector<std::vector<double>> const& GetDistributionSubGrid()   const { return _distributionSubGrid; }   //!< The distribution on the subgrids
    std::vector<double>              const& GetDistributionJointGrid() const { return _distributionJointGrid; } //!< The distribution on the joint grid
    ///@}

  protected:
    Grid                             const& _grid;                  //!< The stored grid reference
    std::vector<std::vector<double>>        _distributionSubGrid;   //!< The array with the distribution values on the subgrid.
    std::vector<double>                     _distributionJointGrid; //!< The array with the distribution values on the joint grid.
  };
}
