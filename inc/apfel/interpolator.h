//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"

#include <vector>
#include <array>
#include <utility>

using std::vector;
using std::array;

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

    Interpolator() = delete;

    /**
     * @brief Interpolator default constructor.
     * @param gr: the x-space grid object over which interpolation takes place
     */
    Interpolator(Grid const& gr);

    /**
     * @name Evaluate functions
     * List of functions that perform the interpolation on the x-space
     * grid.
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
    ///@}

    /**
     * @brief Pure virtual method for the interpolating functions.
     * @param beta: the x-space grid index
     * @param lnx: the value (of the log) of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the interpolation weights
     */
    virtual double Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const = 0;

    /**
     * @brief This purely virtual function computes the lower and
     * upper bounds on which the the sum over interpolants is limited.
     * @param x: the value in x to be interpolated
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the lower and upper bounds of the grid index
     */
    virtual array<int,2> SumBounds(double const& x, SubGrid const& sg) const = 0;

    /**
     * @name Getters
     */
    ///@{
    Grid                   const& GetGrid()                  const { return _grid; }                  //!< The grid
    vector<vector<double>> const& GetDistributionSubGrid()   const { return _distributionSubGrid; }   //!< The distribution on the subgrids
    vector<double>         const& GetDistributionJointGrid() const { return _distributionJointGrid; } //!< The distribution on the joint grid
    ///@}
  protected:
    Grid                   const& _grid;                  //!< The stored grid reference
    vector<vector<double>>        _distributionSubGrid;   //!< The array with the distribution values on the subgrid.
    vector<double>                _distributionJointGrid; //!< The array with the distribution values on the joint grid.
  };
}
