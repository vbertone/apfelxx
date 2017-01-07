//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>
using std::vector;
#include <utility>
using std::pair;

namespace apfel
{
  class Grid;
  class SubGrid;

  /**
   * @brief The Interpolator abstract class
   *
   * This class provides the evaluate interpolation method
   * and requires the implementation of a specialized interpolant algorithm.
   *
   * The current version uses the joint grid object stored allocated by the Grid class.
   */
  class Interpolator
  {
  public:

    Interpolator() = delete;

    /**
     * @brief Interpolator constructor.
     * @param gr the reference to the input grid object
     */
    Interpolator(Grid const& gr);

    /**
     * @brief Evaluate the interpolated function on the joint grid.
     * @param x the requested support value.
     * @return the interpolated result.
     */
    double Evaluate(double const& x) const;

    /**
     * @brief Evaluate the interpolated function on the subgrid with index ig.
     * @param x  the requested support value.
     * @param ig the requested support value.
     * @return the interpolated result.
     */
    double Evaluate(double const& x, int const& ig) const;

    /**
     * @brief Pure virtual method to be defined in the inherited class.
     * @param beta the grid index
     * @param x the value of the required interpolation
     * @param sg SubGrid on which the interpolant is defined
     * @return the interpolation weights.
     */
    virtual double Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const = 0;

    /**
     * @brief Computes the lower and upper bounds on which the the sum over interpolants is limited
     * @param beta the grid index
     * @param x the value of the required interpolation
     * @param sg SubGrid on which the interpolant is defined
     * @return the lower and upper bounds of beta.
     */
    virtual pair<int,int> SumBounds(double const& x, SubGrid const& sg) const = 0;

    // Getters
    Grid                   const& GetGrid()                  const { return _grid; }                  //!< return the grid
    vector<vector<double>> const& GetDistributionSubGrid()   const { return _distributionSubGrid; }   //!< return the distribution on the subgrids
    vector<double>         const& GetDistributionJointGrid() const { return _distributionJointGrid; } //!< return the distribution on the joint grid

  protected:
    Grid                   const& _grid;                  //!< The stored grid reference
    vector<vector<double>>        _distributionSubGrid;   //!< The array with the distribution values on the subgrid.
    vector<double>                _distributionJointGrid; //!< The array with the distribution values on the joint grid.
  };
}
