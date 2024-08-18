//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubledistribution.h"

#include <sstream>

namespace apfel
{
  //_________________________________________________________________________
  DoubleDistribution::DoubleDistribution(Grid                                                const& g1,
                                         Grid                                                const& g2,
                                         std::function<double(double const&, double const&)> const& InDistFunc):
    _g1(g1),
    _g2(g2),
    _li1(LagrangeInterpolator{_g1}),
    _li2(LagrangeInterpolator{_g1})
  {
    const std::vector<double>& jg1 = _g1.GetJointGrid().GetGrid();
    const std::vector<double>& jg2 = _g2.GetJointGrid().GetGrid();

    _dDJointGrid.resize(jg1.size(), jg2.size());
    for (int ix1 = 0; ix1 < (int) jg1.size(); ix1++)
      for (int ix2 = 0; ix2 < (int) jg2.size(); ix2++)
        _dDJointGrid(ix1, ix2) = InDistFunc(std::min(jg1[ix1], 1.), std::min(jg2[ix2], 1.));

    _dDSubGrid.resize(_g1.nGrids());
    for (int ig1 = 0; ig1 < (int) _dDSubGrid.size(); ig1++)
      {
        const std::vector<double>& sg1 = _g1.GetSubGrid(ig1).GetGrid();
        _dDSubGrid[ig1].resize(_g2.nGrids());
        for (int ig2 = 0; ig2 < (int) _dDSubGrid.size(); ig2++)
          {
            const std::vector<double>& sg2 = _g2.GetSubGrid(ig2).GetGrid();
            _dDSubGrid[ig1][ig2].resize(sg1.size(), sg2.size());
            for (int ix1 = 0; ix1 < (int) sg1.size(); ix1++)
              for (int ix2 = 0; ix2 < (int) sg2.size(); ix2++)
                _dDSubGrid[ig1][ig2](ix1, ix2) = InDistFunc(std::min(sg1[ix1], 1.), std::min(sg2[ix2], 1.));
          }
      }
  }

  //_________________________________________________________________________
  DoubleDistribution::DoubleDistribution(DoubleDistribution const& obj):
    _g1(obj.GetFirstGrid()),
    _g2(obj.GetSecondGrid()),
    _li1(LagrangeInterpolator{_g1}),
    _li2(LagrangeInterpolator{_g1}),
    _dDSubGrid(obj.GetDistributionSubGrid()),
    _dDJointGrid(obj.GetDistributionJointGrid())
  {
  }

  //_________________________________________________________________________
  DoubleDistribution::DoubleDistribution(Distribution const& d1, Distribution const& d2):
    _g1(d1.GetGrid()),
    _g2(d2.GetGrid()),
    _li1(LagrangeInterpolator{_g1}),
    _li2(LagrangeInterpolator{_g1})
  {
    const std::vector<double> jg1 = d1.GetDistributionJointGrid();
    const std::vector<double> jg2 = d2.GetDistributionJointGrid();
    _dDJointGrid.resize(jg1.size(), jg2.size());
    for (int ix1 = 0; ix1 < (int) jg1.size(); ix1++)
      for (int ix2 = 0; ix2 < (int) jg2.size(); ix2++)
        _dDJointGrid(ix1, ix2) = jg1[ix1] * jg2[ix2];

    const std::vector<std::vector<double>> sg1 = d1.GetDistributionSubGrid();
    const std::vector<std::vector<double>> sg2 = d2.GetDistributionSubGrid();
    _dDSubGrid.resize(_g1.nGrids());
    for (int ig1 = 0; ig1 < (int) _dDSubGrid.size(); ig1++)
      {
        _dDSubGrid[ig1].resize(_g2.nGrids());
        for (int ig2 = 0; ig2 < (int) _dDSubGrid.size(); ig2++)
          {
            _dDSubGrid[ig1][ig2].resize(sg1.size(), sg2.size());
            for (int ix1 = 0; ix1 < (int) sg1.size(); ix1++)
              for (int ix2 = 0; ix2 < (int) sg2.size(); ix2++)
                _dDSubGrid[ig1][ig2](ix1, ix2) = sg1[ig1][ix1] * sg2[ig2][ix2];
          }
      }
  }

  //_________________________________________________________________________________
  double DoubleDistribution::Evaluate(double const& x1, double const& x2) const
  {
    // Get summation bounds
    const std::array<int, 2> bounds1 = _li1.SumBounds(x1, _g1.GetJointGrid());
    const std::array<int, 2> bounds2 = _li2.SumBounds(x2, _g2.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> intp1(bounds1[1] - bounds1[0]);
    std::vector<double> intp2(bounds2[1] - bounds2[0]);
    for (int beta = 0; beta < bounds1[1] - bounds1[0]; beta++)
      intp1[beta] = _li1.Interpolant(bounds1[0] + beta, x1, _g1.GetJointGrid());
    for (int delta = 0; delta < bounds2[1] - bounds2[0]; delta++)
      intp2[delta] = _li2.Interpolant(bounds2[0] + delta, x2, _g2.GetJointGrid());

    // Interpolate
    double result = 0;
    for (int beta = 0; beta < bounds1[1] - bounds1[0]; beta++)
      for (int delta = 0; delta < bounds2[1] - bounds2[0]; delta++)
        result += intp1[beta] * intp2[delta] * _dDJointGrid(beta + bounds1[0], delta + bounds2[0]);

    return result;
  }

  //_________________________________________________________________________________
  double DoubleDistribution::Derive(double const& x1, double const& x2) const
  {
    // Get summation bounds
    const std::array<int, 2> bounds1 = _li1.SumBounds(x1, _g1.GetJointGrid());
    const std::array<int, 2> bounds2 = _li2.SumBounds(x2, _g2.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> dintp1(bounds1[1] - bounds1[0]);
    std::vector<double> dintp2(bounds2[1] - bounds2[0]);
    for (int beta = 0; beta < bounds1[1] - bounds1[0]; beta++)
      dintp1[beta] = _li1.DerInterpolant(bounds1[0] + beta, x1, _g1.GetJointGrid());
    for (int delta = 0; delta < bounds2[1] - bounds2[0]; delta++)
      dintp2[delta] = _li2.DerInterpolant(bounds2[0] + delta, x2, _g2.GetJointGrid());

    // Interpolate
    double result = 0;
    for (int beta = 0; beta < bounds1[1] - bounds1[0]; beta++)
      for (int delta = 0; delta < bounds2[1] - bounds2[0]; delta++)
        result += dintp1[beta] * dintp2[delta] * _dDJointGrid(beta + bounds1[0], delta + bounds2[0]);

    return result;
  }

  //_________________________________________________________________________________
  double DoubleDistribution::Integrate(double const& a1, double const& b1, double const& a2, double const& b2) const
  {
    // Order integration bounds and adjust the sign if necessary
    const double ao1  = std::min(a1, b1);
    const double bo1  = std::max(a1, b1);
    const double ao2  = std::min(a2, b2);
    const double bo2  = std::max(a2, b2);
    const int    sgn1 = (b1 > a1 ? 1 : -1);
    const int    sgn2 = (b2 > a2 ? 1 : -1);

    // Get summation bounds
    const std::array<int, 2> boundsa1 = _li1.SumBounds(ao1, _g1.GetJointGrid());
    const std::array<int, 2> boundsb1 = _li1.SumBounds(bo1, _g1.GetJointGrid());
    const std::array<int, 2> boundsa2 = _li2.SumBounds(ao2, _g2.GetJointGrid());
    const std::array<int, 2> boundsb2 = _li2.SumBounds(bo2, _g2.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> iintp1(boundsb1[1] - boundsa1[0]);
    std::vector<double> iintp2(boundsb2[1] - boundsa2[0]);
    for (int beta = 0; beta < boundsb1[1] - boundsa1[0]; beta++)
      iintp1[beta] = _li1.IntInterpolant(boundsa1[0] + beta, ao1, bo1, _g1.GetJointGrid());
    for (int delta = 0; delta < boundsb2[1] - boundsa2[0]; delta++)
      iintp2[delta] = _li2.IntInterpolant(boundsa2[0] + delta, ao2, bo2, _g2.GetJointGrid());

    // Interpolate
    double result = 0;
    for (int beta = 0; beta < boundsb1[1] - boundsa1[0]; beta++)
      for (int delta = 0; delta < boundsb2[1] - boundsa2[0]; delta++)
        result += iintp1[beta] * iintp2[delta] * _dDJointGrid(boundsa1[0] + beta, delta + boundsa2[0]);

    return sgn1 * sgn2 * result;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, DoubleDistribution const& in)
  {
    const std::vector<std::vector<matrix<double>>> dd = in.GetDistributionSubGrid();
    os << "DoubleDistribution: " << &in << "\n";
    os << "DoubleDistribution on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) dd.size(); i++)
      for (int j = 0; j < (int) dd[i].size(); j++)
        {
          os << "D[" << i << ", " << j << "]: [";
          for (int alpha = 0; alpha < (int) dd[i][j].size(0); alpha++)
            for (int beta = 0; beta < (int) dd[i][j].size(1); beta++)
              os << "{(" << alpha << ", " << beta << ") : " << dd[i][j](alpha, beta) << "} ";
          os << "]\n";
        }
    os.copyfmt(default_format);
    return os;
  }
}
