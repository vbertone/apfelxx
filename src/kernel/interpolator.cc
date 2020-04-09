//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/interpolator.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  Interpolator::Interpolator(Grid const& gr):
    _grid(gr),
    _distributionJointGrid(_grid.GetJointGrid().GetGrid().size())
  {
    _distributionSubGrid.resize(_grid.nGrids());
    for (int ig = 0; ig < (int) _distributionSubGrid.size(); ig++)
      _distributionSubGrid[ig].resize(_grid.GetSubGrid(ig).GetGrid().size(), 0.);
  }

  //_________________________________________________________________________________
  Interpolator::Interpolator(Grid                             const& gr,
                             std::vector<std::vector<double>> const& distsubgrid,
                             std::vector<double>              const& distjointgrid):
    _grid(gr),
    _distributionSubGrid(distsubgrid),
    _distributionJointGrid(distjointgrid)
  {
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x) const
  {
    const auto bounds = SumBounds(x, _grid.GetJointGrid());
    const double lnx  = log(x);

    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, lnx, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x, int const& ig) const
  {
    const auto bounds = SumBounds(x, _grid.GetSubGrid(ig));
    const double lnx  = log(x);

    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, lnx, _grid.GetSubGrid(ig)) * _distributionSubGrid[ig][beta];

    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Derive(double const& x) const
  {
    const auto bounds = SumBounds(x, _grid.GetJointGrid());
    const double lnx  = log(x);

    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += DerInterpolant(beta, lnx, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    // The factor 1 / x is due to the fact that the original
    // derivative is w.r.t. ln(x) and thus one needs to multiply by
    // dln(x)/dx = 1/x to obtain the derivative w.r.t. x.
    return result / x;
  }

  //_________________________________________________________________________________
  double Interpolator::Integrate(double const& a, double const& b) const
  {
    const auto boundsa = SumBounds(a, _grid.GetJointGrid());
    const auto boundsb = SumBounds(b, _grid.GetJointGrid());
    const double lna  = log(a);
    const double lnb  = log(b);

    double result = 0;
    for (int beta = boundsa[0]; beta < boundsa[1]; beta++)
      result -= IntInterpolant(beta, lna, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    for (int beta = boundsb[0]; beta < boundsb[1]; beta++)
      result += IntInterpolant(beta, lnb, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    return result;
  }
}
