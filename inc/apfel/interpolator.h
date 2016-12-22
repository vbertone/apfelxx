/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@cern.ch
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */
#pragma once

#include <vector>
using std::vector;

namespace apfel {

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
    /**
     * @brief Interpolator constructor.
     * @param gr the reference to the input grid object
     */
    Interpolator(Grid const& gr);

    /**
     * @brief Evaluate the interpolated function by the grid.
     * @param x the requested support value.
     * @return the interpolated result.
     */
    double Evaluate(double const& x) const;

    /**
     * @brief Pure virtual method to be defined in the inherited class.
     * @param beta the grid index
     * @param x the value of the required interpolation
     * @param sg the subgrid object.
     * @return the interpolation weights.
     */
    virtual double Interpolant(int const& beta, double const& x, SubGrid const& sg) const = 0;

  protected:
    Grid const& _grid;            //!< The stored grid reference
    vector<double> _distribution; //!< The array with the distribution values for the joint grid.
  };

  /**
   * @brief The LagrangeInterpolator class.
   *
   * A specialization example of the Interpolator
   * class using the lagrange interpolation.
   */
  class LagrangeInterpolator: public Interpolator
  {
  public:

    /**
     * @see Interpolator::Interpolator
     */
    LagrangeInterpolator(Grid const& gr);

    /**
     * @see Interpolator::Interpolant
     */
    double Interpolant(int const& beta, const double &x, SubGrid const& sg) const;
  };
}
