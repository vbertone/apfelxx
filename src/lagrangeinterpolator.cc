//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <cmath>
#include <iostream>

#include "apfel/lagrangeinterpolator.h"
#include "apfel/grid.h"
#include "apfel/subgrid.h"
#include "apfel/tools.h"

using namespace std;

namespace apfel {

  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(Grid const& gr):
    Interpolator{gr}
  {
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::Interpolant(int const& beta, double const& lnx, SubGrid const& sg) const
  {
    // Compute interpolant
    auto const& lxsg = sg.GetLogGrid();

    // Return immediately 1 is "x" coincides with "xg[beta]".
    if (fabs(lnx - lxsg[beta]) < eps12) return 1;

    // Return 0 if "x" is outside the range in which the interpolant is different from zero.
    // Ideally this functions should never be called if "beta" and "x" are such that "Interpolant"
    // is identically zero. Use "SumBounds" to know where "beta" should run over given "x".
    auto const id = sg.InterDegree();
    int bound = beta - id;
    if (id > beta) bound = 0;
    if (lnx < lxsg[bound] || lnx >= lxsg[beta+1]) return 0;

    // Initialize interpolant
    double w_int = 1;

    // Find the the neighbors of "x" on the grid
    int j;
    for (j = 0; j <= beta-bound; j++)
      if (lnx >= lxsg[beta-j] && lnx < lxsg[beta-j+1])
	break;

    // Compute the interpolant
    for (auto delta = 0; delta <= id; delta++)
      if(delta != j)
	w_int *= ( lnx - lxsg[beta-j+delta] )  / ( lxsg[beta] - lxsg[beta-j+delta] );

    return w_int;   
  }

  //_________________________________________________________________________________
  pair<int,int> LagrangeInterpolator::SumBounds(double const& x, SubGrid const& sg) const
  {
    auto const& xsg = sg.GetGrid();

    pair<int,int> bounds (0,0);
    if (x < xsg[0] - eps12 || x > xsg[sg.nx()] + eps12)
      return bounds;

    const auto low = lower_bound(xsg.begin()+1, xsg.end()-sg.InterDegree()-1, x) - xsg.begin();
    bounds = {low, low};

    if (fabs(x - xsg[low]) <= eps12)
          bounds.second += 1;
    else
      {
        bounds.first -= 1;
        bounds.second += sg.InterDegree();
      }

    return bounds;
  }

}
