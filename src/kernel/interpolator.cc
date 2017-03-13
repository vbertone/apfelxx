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
    auto const bounds = SumBounds(x, _grid.GetJointGrid());
    auto const lnx    = log(x);

    double result = 0;
    for (auto beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, lnx, _grid.GetJointGrid()) * _distributionJointGrid[beta];
    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x, int const& ig) const
  {
    auto const bounds = SumBounds(x, _grid.GetSubGrid(ig));
    auto const lnx    = log(x);

    double result = 0;
    for (auto beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, lnx, _grid.GetSubGrid(ig)) * _distributionSubGrid[ig][beta];
    return result;
  }

}
