//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/interpolator.h"

namespace apfel
{
  /**
   * @brief The LagrangeInterpolator class is a specialization of the
   * Interpolator class using the lagrange interpolation procedure.
   */
  class LagrangeInterpolator: public Interpolator
  {
  public:
    /**
     * @brief The LagrangeInterpolator constructor.
     * @param gr: the x-space grid object over which interpolation takes place
     * @see Interpolator::Interpolator
     */
    LagrangeInterpolator(Grid const& gr);

    /**
     * @brief The LagrangeInterpolator constructor.
     * @param gr: the x-space grid object over which interpolation takes place
     * @param distsubgrid: the vector of subgrids
     * @param distjointgrid: the joint subgrid
     * @see Interpolator::Interpolator
     * @note distjointgrid and distsubgrid are assumed to match the
     * structure of the grid gr.
     */
    LagrangeInterpolator(Grid                             const& gr,
                         std::vector<std::vector<double>> const& distsubgrid,
                         std::vector<double>              const& distjointgrid);

    /**
     * @brief This function defines the interpolating function used by
     * the mothe class Interpolator to perform the actual
     * interpolation.
     * @param beta: the x-space grid index
     * @param lnx: the value (of the log) of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the interpolation weights
     * @see Interpolator::Interpolant
     */
    double Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const;

    /**
     * @brief This function computes the lower and upper bounds on
     * which the the sum over interpolants is limited.
     * @param x: the value in x to be interpolated
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the lower and upper bounds of the grid index
     * @see Interpolator::SumBounds
     */
    std::array<int,2> SumBounds(double const& x, SubGrid const& sg) const;
  };
}
