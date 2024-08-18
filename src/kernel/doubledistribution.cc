//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubledistribution.h"

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
    const std::array<int, 2> bounds1 = _li1.SumBounds(x1, _g1.GetJointGrid());
    const std::array<int, 2> bounds2 = _li2.SumBounds(x2, _g2.GetJointGrid());
    double result = 0;
    for (int beta = bounds1[0]; beta < bounds1[1]; beta++)
      {
        const double intp1 = _li1.Interpolant(beta, x1, _g1.GetJointGrid());
        for (int delta = bounds2[0]; delta < bounds2[1]; delta++)
          {
            const double intp2 = _li2.Interpolant(delta, x2, _g2.GetJointGrid());
            result += intp1 * intp2 * _dDJointGrid(beta, delta);
          }
      }
    return result;
  }

  //_________________________________________________________________________________
  double DoubleDistribution::Evaluate(double const& x1, double const& x2, int const& ig1, int const& ig2) const
  {
    const std::array<int, 2> bounds1 = _li1.SumBounds(x1, _g1.GetJointGrid());
    const std::array<int, 2> bounds2 = _li2.SumBounds(x2, _g2.GetJointGrid());
    double result = 0;
    for (int beta = bounds1[0]; beta < bounds1[1]; beta++)
      {
        const double intp1 = _li1.Interpolant(beta, x1, _g1.GetJointGrid());
        for (int delta = bounds2[0]; delta < bounds2[1]; delta++)
          {
            const double intp2 = _li2.Interpolant(delta, x2, _g2.GetJointGrid());
            result += intp1 * intp2 * _dDSubGrid[ig1][ig2](beta, delta);
          }
      }
    return result;
  }
}
