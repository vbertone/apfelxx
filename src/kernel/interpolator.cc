//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/interpolator.h"
#include "apfel/messages.h"

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
    const std::array<int, 2> bounds = SumBounds(x, _grid.GetJointGrid());
    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, x, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Evaluate(double const& x, int const& ig) const
  {
    const std::array<int, 2> bounds = SumBounds(x, _grid.GetSubGrid(ig));
    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += Interpolant(beta, x, _grid.GetSubGrid(ig)) * _distributionSubGrid[ig][beta];

    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Derive(double const& x) const
  {
    const std::array<int, 2> bounds = SumBounds(x, _grid.GetJointGrid());
    double result = 0;
    for (int beta = bounds[0]; beta < bounds[1]; beta++)
      result += DerInterpolant(beta, x, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    return result;
  }

  //_________________________________________________________________________________
  double Interpolator::Integrate(double const& a, double const& b) const
  {
    // Order integration bounds and adjust the sign if necessary
    const double ao  = std::min(a, b);
    const double bo  = std::max(a, b);
    const int    sgn = (b > a ? 1 : -1);

    // Get summation bounds
    const std::array<int, 2> boundsa = SumBounds(ao, _grid.GetJointGrid());
    const std::array<int, 2> boundsb = SumBounds(bo, _grid.GetJointGrid());

    // Interpolate
    double result = 0;
    for (int beta = boundsa[0]; beta < boundsb[1]; beta++)
      result += IntInterpolant(beta, ao, bo, _grid.GetJointGrid()) * _distributionJointGrid[beta];

    return sgn * result;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, Interpolator const& in)
  {
    const std::vector<double> sb = in.GetDistributionSubGrid()[0];
    os << "Interpolator: " << &in << "\n";
    os << "Distribution on the first SubGrid:\n[";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(1);
    for (int i = 0; i < (int) sb.size(); i++)
      os << sb[i] << " ";
    os << "\b]";
    os.copyfmt(default_format);
    return os;
  }
}
