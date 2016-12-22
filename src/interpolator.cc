/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@cern.ch
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */

#include "apfel/interpolator.h"
#include "apfel/grid.h"
#include "apfel/subgrid.h"
#include <cmath>
#include <iostream>

namespace apfel {

  //_________________________________________________________________________________
  Interpolator::Interpolator(Grid const& gr):
    _grid(gr)
  {
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x) const
  {
    double result = 0;
    for (auto beta = 0; beta < _grid.GetJointGrid().nx(); beta++)
      result += Interpolant(beta, x, _grid.GetJointGrid())*_distribution[beta];
    return result;
  }

  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(const Grid &gr):
    Interpolator{gr}
  {
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::Interpolant(const int &beta, const double &x, const SubGrid &sg) const
  {    
    // Compute interpolant
    const auto xsg = sg.GetGrid();
    const auto lxsg = sg.GetLogGrid();
    int bound = beta - sg.InterDegree();
    if(sg.InterDegree() > beta) bound = 0;
    if(x < xsg[bound] || x >= xsg[beta+1]) return 0;

    double w_int = 1;

    // If the grid is external ...
    if(sg.IsExternal())
      {
        int j;
        for(j = 0; j <= beta-bound; j++)
          if(x >= xsg[beta-j] && x < xsg[beta-j+1]) break;

        for(auto delta = 0; delta <= sg.InterDegree(); delta++)
          if(delta != j)
            w_int *= (log(x) - lxsg[beta-j+delta])  / ( lxsg[beta] - lxsg[beta-j+delta] );
      }
    else  // If the grid is internal ...
      {
        double fact = (log(x) -lxsg[beta] ) / sg.Step();
        int j;
        for(j=0; j<=beta-bound; j++)
          if(x >= xsg[beta-j] && x < xsg[beta-j+1])
            break;

        for(auto delta = 0; delta <= sg.InterDegree(); delta++)
          if(delta != j)
            w_int *= ( fact / ( j - delta ) + 1 );
      }

    return w_int;
  }
}
