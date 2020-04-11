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
     * the mother class Interpolator to perform the actual
     * interpolation using polynomials in log(x).
     * @param beta: the x-space grid index
     * @param lnx: the value (of the log) of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the interpolation weights
     * @see Interpolator::InterpolantLog
     */
    double InterpolantLog(int const& beta, double const& lnx, SubGrid const& sg) const;

    /**
     * @brief This function defines the interpolating function used by
     * the mother class Interpolator to perform the
     * interpolation.
     * @param beta: the x-space grid index
     * @param x: the value of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the interpolation weights
     * @see Interpolator::Interpolant
     */
    double Interpolant(int const& beta, double const& x, SubGrid const& sg) const;

    /**
     * @brief This function defines the derivative of the
     * interpolating function used by the mother class Interpolator to
     * perform the actual interpolation.
     * @param beta: the x-space grid index
     * @param x: the value of the interpolation point
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the derivative of the interpolation weights
     * @see Interpolator::DerInterpolant
     */
    double DerInterpolant(int const& beta, double const& x, SubGrid const& sg) const;

    /**
     * @brief This function defines the integral of the interpolating
     * function used by the mother class Interpolator to perform the
     * actual interpolation.
     * @param beta: the x-space grid index
     * @param a: the value of the lower integration bound
     * @param b: the value of the upper integration bound
     * @param sg: the SubGrid over which the interpolant is defined
     * @return the integral of the interpolation weights
     * @see Interpolator::IntInterpolant
     */
    double IntInterpolant(int const& beta, double const& a, double const& b, SubGrid const& sg) const;

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
