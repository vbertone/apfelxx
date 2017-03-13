//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/lagrangeinterpolator.h"
#include "apfel/tools.h"

#include <cmath>
#include <iostream>

namespace apfel {

  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(Grid const& gr):
    Interpolator{gr}
  {
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const
  {
    // Get the logarithmic grid
    auto const& lxsg = sg.GetLogGrid();

    // Return immediately 1 if "x" coincides with "xg[beta]".
    if (fabs(lnx - lxsg[beta]) < eps12) return 1;

    // Return 0 if "x" is outside the range in which the interpolant is different from zero.
    // Ideally this functions should never be called if "beta" and "x" are such that "Interpolant"
    // is identically zero. Use "SumBounds" to know where "beta" should run over given "x".
    auto const id = sg.InterDegree();
    int bound = beta - id;
    if (id > beta) bound = 0;
    if (lnx < lxsg[bound] || lnx >= lxsg[beta+1]) return 0;

    // Find the the neighbors of "x" on the grid
    int j;
    for (j = 0; j <= beta-bound; j++)
      if (lnx >= lxsg[beta-j])
	break;

    // Compute the interpolant
    double w_int = 1;
    for (auto delta = beta-j; delta <= beta-j+id; delta++)
      if (delta != beta)
        w_int *= ( lnx - lxsg[delta] ) / ( lxsg[beta] - lxsg[delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  array<int, 2> LagrangeInterpolator::SumBounds(double const& x, SubGrid const& sg) const
  {
    auto const& xsg = sg.GetGrid();

    array<int,2> bounds = {0, 0};
    if (x < xsg[0] - eps12 || x > xsg[sg.nx()] + eps12)
      return bounds;

    const auto low = lower_bound(xsg.begin()+1, xsg.end()-sg.InterDegree()-1, x) - xsg.begin();
    bounds[0] = bounds[1] = low;

    if (fabs(x - xsg[low]) <= eps12)
          bounds[1] += 1;
    else
      {
        bounds[0] -= 1;
        bounds[1] += sg.InterDegree();
      }

    return bounds;
  }

}
