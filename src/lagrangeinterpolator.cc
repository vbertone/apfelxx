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
  LagrangeInterpolator::LagrangeInterpolator(const Grid &gr):
    Interpolator{gr}
  {
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::Interpolant(const int &beta, double const& x, SubGrid const& sg) const
  {
    // Compute interpolant
    const auto xsg  = sg.GetGrid();

    // Return immediately 1 is "x" coincides with "xg[beta]".
    if (abs(x - xsg[beta]) < eps12) return 1;

    // Return 0 if "x" is outside the range in which the interpolant is different from zero.
    // Ideally this functions should never be called if "beta" and "x" are such that "Interpolant"
    // is identically zero. Use "SumBounds" to know where "beta" should run over given "x".
    const auto id   = sg.InterDegree();
    int bound = beta - id;
    if (id > beta) bound = 0;
    if (x < xsg[bound] || x >= xsg[beta+1]) return 0;

    // Initialize interpolant
    double w_int = 1;

    // Find the the neighbors of "x" on the grid
    int j;
    for (j = 0; j <= beta-bound; j++)
      if (x >= xsg[beta-j] && x < xsg[beta-j+1])
	break;

    const auto lxsg = sg.GetLogGrid();
    double lnx = log(x);

    // If the grid is external ...
    if(sg.IsExternal())
      {
        for (auto delta = 0; delta <= id; delta++)
          if(delta != j)
            w_int *= ( lnx - lxsg[beta-j+delta] )  / ( lxsg[beta] - lxsg[beta-j+delta] );
      }
    // If the grid is internal ...
    else
      {
        double fact = ( lnx - lxsg[beta] ) / sg.Step();
        for (auto delta = 0; delta <= id; delta++)
          if(delta != j)
            w_int *= ( fact / ( j - delta ) + 1 );
      }
    return w_int;
  }

  //_________________________________________________________________________________
  pair<int,int> LagrangeInterpolator::SumBounds(double const& x, SubGrid const& sg) const
  {
    const auto xsg = sg.GetGrid();    
    const auto n   = sg.nx();
    const auto id = sg.InterDegree();

    pair<int,int> bounds;
    for (auto beta = 1; beta <= n; beta++)
      {
        // If "x" coincides with "xsg[beta]" within a certain accuracy, the interpolation
        // function is delta e thus it is enough to sum only over one value of "beta"...
        if (fabs(x - xsg[beta]) < eps12)
          {
            bounds.first  = beta;
            bounds.second = beta + 1;
          }
        // ... othewise if the range extends from the lower neighbor of "x" on the grid
        // to "id" more nodes on the right of "x".
        else if (xsg[beta] > x)
          {
            bounds.first  = beta - 1;
            bounds.second = beta + id;
            break;
          }
      }

    return bounds;
  }

}
