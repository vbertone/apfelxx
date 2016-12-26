//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

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
    auto const bounds = SumBounds(x, _grid.GetJointGrid());
    double result = 0;
    for (auto beta = bounds.first; beta < bounds.second; beta++)
      result += Interpolant(beta, x, _grid.GetJointGrid()) * _distributionJointGrid[beta];
    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x, int const& ig) const
  {
    auto const bounds = SumBounds(x, _grid.GetSubGrid(ig));
    double result = 0;
    for (auto beta = bounds.first; beta < bounds.second; beta++)
      result += Interpolant(beta, x, _grid.GetSubGrid(ig)) * _distributionSubGrid[ig][beta];
    return result;
  }

}
