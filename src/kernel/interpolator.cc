//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/interpolator.h"

#include <cmath>

namespace apfel {
  //_________________________________________________________________________________
  Interpolator::Interpolator(Grid const& gr):
    _grid(gr)
  {
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x) const
  {
    const auto   bounds = SumBounds(x, _grid.GetJointGrid());
    const double lnx    = log(x);

    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, lnx, _grid.GetJointGrid()) * _distributionJointGrid[beta];
    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x, int const& ig) const
  {
    const auto   bounds = SumBounds(x, _grid.GetSubGrid(ig));
    const double lnx    = log(x);

    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, lnx, _grid.GetSubGrid(ig)) * _distributionSubGrid[ig][beta];
    return result;
  }
}
