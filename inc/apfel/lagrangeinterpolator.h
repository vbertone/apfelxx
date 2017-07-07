//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/interpolator.h"

namespace apfel
{
  /**
   * @brief The LagrangeInterpolator class.
   *
   * A specialization example of the Interpolator class using the
   * lagrange interpolation.
   */
  class LagrangeInterpolator: public Interpolator
  {
  public:

    LagrangeInterpolator() = delete;

    /**
     * @see Interpolator::Interpolator
     */
    LagrangeInterpolator(Grid const& gr);

    /**
     * @see Interpolator::Interpolant
     */
    double Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const;

    /**
     * @see Interpolator::SumBounds
     */
    array<int,2> SumBounds(double const& x, SubGrid const& sg) const;
  };
}
