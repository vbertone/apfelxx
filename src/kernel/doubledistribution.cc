//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubledistribution.h"
#include "apfel/messages.h"
#include "apfel/constants.h"

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
    _li2(LagrangeInterpolator{_g2})
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
    _li2(LagrangeInterpolator{_g2}),
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
  DoubleDistribution::DoubleDistribution(Grid                                     const& g1,
                                         Grid                                     const& g2,
                                         std::vector<std::vector<matrix<double>>> const& distsubgrid,
                                         matrix<double>                           const& distjointgrid):
    _g1(g1),
    _g2(g2),
    _li1(LagrangeInterpolator{_g1}),
    _li2(LagrangeInterpolator{_g2}),
    _dDSubGrid(distsubgrid),
    _dDJointGrid(distjointgrid)
  {
  }

  //_________________________________________________________________________________
  DoubleDistribution::DoubleDistribution(DoubleObject<Distribution> const& DObj):
    _g1(DObj.GetTerms()[0].object1.GetGrid()),
    _g2(DObj.GetTerms()[0].object2.GetGrid()),
    _li1(LagrangeInterpolator{_g1}),
    _li2(LagrangeInterpolator{_g2})
  {
    // Number on nodes on the joint grids
    const int nxj1 = _g1.GetJointGrid().GetGrid().size();
    const int nxj2 = _g2.GetJointGrid().GetGrid().size();

    // Number of subgrids
    const int ng1 = _g1.nGrids();
    const int ng2 = _g2.nGrids();

    // Resize matrices
    _dDJointGrid.resize(nxj1, nxj2);
    _dDSubGrid.resize(ng1, std::vector<matrix<double>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        _dDSubGrid[ig1][ig2].resize(_g1.GetSubGrid(ig1).GetGrid().size(), _g2.GetSubGrid(ig2).GetGrid().size());

    // Fill in matrices
    for (auto const& term : DObj.GetTerms())
      {
        // Get coefficient and objects
        const double coef = term.coefficient;
        const std::vector<double> Djg1 = term.object1.GetDistributionJointGrid();
        const std::vector<double> Djg2 = term.object2.GetDistributionJointGrid();
        const std::vector<std::vector<double>> Dsg1 = term.object1.GetDistributionSubGrid();
        const std::vector<std::vector<double>> Dsg2 = term.object2.GetDistributionSubGrid();

        // Fill in joint matrix
        for (int i1 = 0; i1 < nxj1; i1++)
          for (int i2 = 0; i2 < nxj2; i2++)
            _dDJointGrid(i1, i2) += coef * Djg1[i1] * Djg2[i2];

        // Fill in subgrids
        for (int ig1 = 0; ig1 < ng1; ig1++)
          for (int ig2 = 0; ig2 < ng2; ig2++)
            for (int i1 = 0; i1 < (int) _dDSubGrid[ig1][ig2].size(0); i1++)
              for (int i2 = 0; i2 < (int) _dDSubGrid[ig1][ig2].size(1); i2++)
                _dDSubGrid[ig1][ig2](i1, i2) += coef * Dsg1[ig1][i1] * Dsg2[ig2][i2];
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

  //_________________________________________________________________________
  Distribution DoubleDistribution::Evaluate1(double const& x1) const
  {
    // Get summation bounds on the joint grid along the first
    // direction.
    const std::array<int, 2> bounds1j = _li1.SumBounds(x1, _g1.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> intp1j(bounds1j[1] - bounds1j[0]);
    for (int beta = 0; beta < bounds1j[1] - bounds1j[0]; beta++)
      intp1j[beta] = _li1.Interpolant(bounds1j[0] + beta, x1, _g1.GetJointGrid());

    // Get number of nodes of the joint grid along the second
    // direction.
    const int nxj2 = _g2.GetJointGrid().GetGrid().size();

    // Interpolate along the first direction of the joint grid
    std::vector<double> jg(nxj2, 0.);
    for (int delta = 0; delta < nxj2; delta++)
      for (int beta = 0; beta < bounds1j[1] - bounds1j[0]; beta++)
        jg[delta] += intp1j[beta] * _dDJointGrid(beta + bounds1j[0], delta);

    // Now take care of the subgrids. ig1 has to correspond to the
    // subgrid index along the first direction such that "x1" is on
    // the denser grid possible.
    int ig1;
    for (ig1 = 0; ig1 < _g1.nGrids(); ig1++)
      if(_g1.GetSubGrid(ig1).xMin() > x1)
        break;
    ig1--;

    // Get summation bounds on the ig1-th subgrid along the first
    // direction.
    const std::array<int, 2> bounds1s = _li1.SumBounds(x1, _g1.GetSubGrid(ig1));

    // Accumulate interpolating functions
    std::vector<double> intp1s(bounds1s[1] - bounds1s[0]);
    for (int beta = 0; beta < bounds1s[1] - bounds1s[0]; beta++)
      intp1s[beta] = _li1.Interpolant(bounds1s[0] + beta, x1, _g1.GetSubGrid(ig1));

    // Run over the subgrids along the second direction
    const int ng2 = _g2.nGrids();
    std::vector<std::vector<double>> sg(ng2);
    for (int ig2 = 0; ig2 < ng2; ig2++)
      {
        // Get number of nodes of the ig2-th subgrid along the second
        // direction.
        const int nx2s = _g2.GetSubGrid(ig2).GetGrid().size();

        // Interpolate along the first direction of the ig2-th subgrid
        sg[ig2].resize(nx2s, 0.);
        for (int delta = 0; delta < nx2s; delta++)
          for (int beta = 0; beta < bounds1s[1] - bounds1s[0]; beta++)
            sg[ig2][delta] += intp1s[beta] * _dDSubGrid[ig1][ig2](beta + bounds1j[0], delta);
      }
    return Distribution{_g2, sg, jg};
  }

  //_________________________________________________________________________
  Distribution DoubleDistribution::Evaluate2(double const& x2) const
  {
    // Get summation bounds on the joint grid along the second
    // direction.
    const std::array<int, 2> bounds2j = _li2.SumBounds(x2, _g2.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> intp2j(bounds2j[1] - bounds2j[0]);
    for (int delta = 0; delta < bounds2j[1] - bounds2j[0]; delta++)
      intp2j[delta] = _li2.Interpolant(bounds2j[0] + delta, x2, _g2.GetJointGrid());

    // Get number of nodes of the joint grid along the first
    // direction.
    const int nxj1 = _g1.GetJointGrid().GetGrid().size();

    // Interpolate along the first direction of the joint grid
    std::vector<double> jg(nxj1, 0.);
    for (int beta = 0; beta < nxj1; beta++)
      for (int delta = 0; delta < bounds2j[1] - bounds2j[0]; delta++)
        jg[beta] += intp2j[delta] * _dDJointGrid(beta, delta + bounds2j[0]);

    // Now take care of the subgrids. ig2 has to correspond to the
    // subgrid index along the second direction such that "x2" is on
    // the denser grid possible.
    int ig2;
    for (ig2 = 0; ig2 < _g2.nGrids(); ig2++)
      if(_g2.GetSubGrid(ig2).xMin() > x2)
        break;
    ig2--;

    // Get summation bounds on the igw-th subgrid along the second
    // direction.
    const std::array<int, 2> bounds2s = _li2.SumBounds(x2, _g2.GetSubGrid(ig2));

    // Accumulate interpolating functions
    std::vector<double> intp2s(bounds2s[1] - bounds2s[0]);
    for (int delta = 0; delta < bounds2s[1] - bounds2s[0]; delta++)
      intp2s[delta] = _li2.Interpolant(bounds2s[0] + delta, x2, _g2.GetSubGrid(ig2));

    // Run over the subgrids along the first direction
    const int ng1 = _g1.nGrids();
    std::vector<std::vector<double>> sg(ng1);
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Get number of nodes of the ig1-th subgrid along the first
        // direction.
        const int nx1s = _g1.GetSubGrid(ig1).GetGrid().size();

        // Interpolate along the second direction of the ig1-th subgrid
        sg[ig1].resize(nx1s, 0.);
        for (int beta = 0; beta < nx1s; beta++)
          for (int delta = 0; delta < bounds2s[1] - bounds2s[0]; delta++)
            sg[ig1][beta] += intp2s[delta] * _dDSubGrid[ig1][ig2](beta, delta + bounds2j[0]);
      }
    return Distribution{_g1, sg, jg};
  }

  //_________________________________________________________________________
  Distribution DoubleDistribution::Derive1(double const& x1) const
  {
    // Get summation bounds on the joint grid along the first
    // direction.
    const std::array<int, 2> bounds1j = _li1.SumBounds(x1, _g1.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> dintp1j(bounds1j[1] - bounds1j[0]);
    for (int beta = 0; beta < bounds1j[1] - bounds1j[0]; beta++)
      dintp1j[beta] = _li1.DerInterpolant(bounds1j[0] + beta, x1, _g1.GetJointGrid());

    // Get number of nodes of the joint grid along the second
    // direction.
    const int nxj2 = _g2.GetJointGrid().GetGrid().size();

    // Interpolate along the first direction of the joint grid
    std::vector<double> jg(nxj2, 0.);
    for (int delta = 0; delta < nxj2; delta++)
      for (int beta = 0; beta < bounds1j[1] - bounds1j[0]; beta++)
        jg[delta] += dintp1j[beta] * _dDJointGrid(beta + bounds1j[0], delta);

    // Now take care of the subgrids. ig1 has to correspond to the
    // subgrid index along the first direction such that "x1" is on
    // the denser grid possible.
    int ig1;
    for (ig1 = 0; ig1 < _g1.nGrids(); ig1++)
      if(_g1.GetSubGrid(ig1).xMin() > x1)
        break;
    ig1--;

    // Get summation bounds on the ig1-th subgrid along the first
    // direction.
    const std::array<int, 2> bounds1s = _li1.SumBounds(x1, _g1.GetSubGrid(ig1));

    // Accumulate interpolating functions
    std::vector<double> dintp1s(bounds1s[1] - bounds1s[0]);
    for (int beta = 0; beta < bounds1s[1] - bounds1s[0]; beta++)
      dintp1s[beta] = _li1.DerInterpolant(bounds1s[0] + beta, x1, _g1.GetSubGrid(ig1));

    // Run over the subgrids along the second direction
    const int ng2 = _g2.nGrids();
    std::vector<std::vector<double>> sg(ng2);
    for (int ig2 = 0; ig2 < ng2; ig2++)
      {
        // Get number of nodes of the ig2-th subgrid along the second
        // direction.
        const int nx2s = _g2.GetSubGrid(ig2).GetGrid().size();

        // Interpolate along the first direction of the ig2-th subgrid
        sg[ig2].resize(nx2s, 0.);
        for (int delta = 0; delta < nx2s; delta++)
          for (int beta = 0; beta < bounds1s[1] - bounds1s[0]; beta++)
            sg[ig2][delta] += dintp1s[beta] * _dDSubGrid[ig1][ig2](beta + bounds1j[0], delta);
      }
    return Distribution{_g2, sg, jg};
  }

  //_________________________________________________________________________
  Distribution DoubleDistribution::Derive2(double const& x2) const
  {
    // Get summation bounds on the joint grid along the second
    // direction.
    const std::array<int, 2> bounds2j = _li2.SumBounds(x2, _g2.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> dintp2j(bounds2j[1] - bounds2j[0]);
    for (int delta = 0; delta < bounds2j[1] - bounds2j[0]; delta++)
      dintp2j[delta] = _li2.DerInterpolant(bounds2j[0] + delta, x2, _g2.GetJointGrid());

    // Get number of nodes of the joint grid along the first
    // direction.
    const int nxj1 = _g1.GetJointGrid().GetGrid().size();

    // Interpolate along the first direction of the joint grid
    std::vector<double> jg(nxj1, 0.);
    for (int beta = 0; beta < nxj1; beta++)
      for (int delta = 0; delta < bounds2j[1] - bounds2j[0]; delta++)
        jg[beta] += dintp2j[delta] * _dDJointGrid(beta, delta + bounds2j[0]);

    // Now take care of the subgrids. ig2 has to correspond to the
    // subgrid index along the second direction such that "x2" is on
    // the denser grid possible.
    int ig2;
    for (ig2 = 0; ig2 < _g2.nGrids(); ig2++)
      if(_g2.GetSubGrid(ig2).xMin() > x2)
        break;
    ig2--;

    // Get summation bounds on the igw-th subgrid along the second
    // direction.
    const std::array<int, 2> bounds2s = _li2.SumBounds(x2, _g2.GetSubGrid(ig2));

    // Accumulate interpolating functions
    std::vector<double> dintp2s(bounds2s[1] - bounds2s[0]);
    for (int delta = 0; delta < bounds2s[1] - bounds2s[0]; delta++)
      dintp2s[delta] = _li2.DerInterpolant(bounds2s[0] + delta, x2, _g2.GetSubGrid(ig2));

    // Run over the subgrids along the first direction
    const int ng1 = _g1.nGrids();
    std::vector<std::vector<double>> sg(ng1);
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Get number of nodes of the ig1-th subgrid along the first
        // direction.
        const int nx1s = _g1.GetSubGrid(ig1).GetGrid().size();

        // Interpolate along the second direction of the ig1-th subgrid
        sg[ig1].resize(nx1s, 0.);
        for (int beta = 0; beta < nx1s; beta++)
          for (int delta = 0; delta < bounds2s[1] - bounds2s[0]; delta++)
            sg[ig1][beta] += dintp2s[delta] * _dDSubGrid[ig1][ig2](beta, delta + bounds2j[0]);
      }
    return Distribution{_g1, sg, jg};
  }

  //_________________________________________________________________________
  Distribution DoubleDistribution::Integrate1(double const& a1, double const& b1) const
  {
    // Order integration bounds and adjust the sign if necessary
    const double ao1  = std::min(a1, b1);
    const double bo1  = std::max(a1, b1);
    const int    sgn1 = (b1 > a1 ? 1 : -1);

    // Get summation bounds on the joint grid along the first
    // direction.
    const std::array<int, 2> boundsa1j = _li1.SumBounds(ao1, _g1.GetJointGrid());
    const std::array<int, 2> boundsb1j = _li1.SumBounds(bo1, _g1.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> iintp1j(boundsb1j[1] - boundsa1j[0]);
    for (int beta = 0; beta < boundsb1j[1] - boundsa1j[0]; beta++)
      iintp1j[beta] = _li1.IntInterpolant(boundsa1j[0] + beta, ao1, bo1, _g1.GetJointGrid());

    // Get number of nodes of the joint grid along the second
    // direction.
    const int nxj2 = _g2.GetJointGrid().GetGrid().size();

    // Interpolate along the first direction of the joint grid
    std::vector<double> jg(nxj2, 0.);
    for (int delta = 0; delta < nxj2; delta++)
      for (int beta = 0; beta < boundsb1j[1] - boundsa1j[0]; beta++)
        jg[delta] += sgn1 * iintp1j[beta] * _dDJointGrid(beta + boundsa1j[0], delta);

    // Now take care of the subgrids. ig1 has to correspond to the
    // subgrid index along the first direction such that the lower
    // bound "ao1" is on the denser grid possible.
    int ig1;
    for (ig1 = 0; ig1 < _g1.nGrids(); ig1++)
      if(_g1.GetSubGrid(ig1).xMin() > ao1)
        break;
    ig1--;

    // Get summation bounds on the ig1-th subgrid along the first
    // direction.
    const std::array<int, 2> boundsa1s = _li1.SumBounds(ao1, _g1.GetSubGrid(ig1));
    const std::array<int, 2> boundsb1s = _li1.SumBounds(bo1, _g1.GetSubGrid(ig1));

    // Accumulate interpolating functions
    std::vector<double> iintp1s(boundsb1s[1] - boundsa1s[0]);
    for (int beta = 0; beta < boundsb1s[1] - boundsa1s[0]; beta++)
      iintp1s[beta] = _li1.IntInterpolant(boundsa1s[0] + beta, ao1, bo1, _g1.GetSubGrid(ig1));

    // Run over the subgrids along the second direction
    const int ng2 = _g2.nGrids();
    std::vector<std::vector<double>> sg(ng2);
    for (int ig2 = 0; ig2 < ng2; ig2++)
      {
        // Get number of nodes of the ig2-th subgrid along the second
        // direction.
        const int nx2s = _g2.GetSubGrid(ig2).GetGrid().size();

        // Interpolate along the first direction of the ig2-th subgrid
        sg[ig2].resize(nx2s, 0.);
        for (int delta = 0; delta < nx2s; delta++)
          for (int beta = 0; beta < boundsb1s[1] - boundsa1s[0]; beta++)
            sg[ig2][delta] += sgn1 * iintp1s[beta] * _dDSubGrid[ig1][ig2](beta + boundsa1j[0], delta);
      }
    return Distribution{_g2, sg, jg};
  }

  //_________________________________________________________________________
  Distribution DoubleDistribution::Integrate2(double const& a2, double const& b2) const
  {
    // Order integration bounds and adjust the sign if necessary
    const double ao2  = std::min(a2, b2);
    const double bo2  = std::max(a2, b2);
    const int    sgn2 = (b2 > a2 ? 1 : -1);

    // Get summation bounds on the joint grid along the second
    // direction.
    const std::array<int, 2> boundsa2j = _li2.SumBounds(ao2, _g2.GetJointGrid());
    const std::array<int, 2> boundsb2j = _li2.SumBounds(bo2, _g2.GetJointGrid());

    // Accumulate interpolating functions
    std::vector<double> iintp2j(boundsb2j[1] - boundsa2j[0]);
    for (int delta = 0; delta < boundsb2j[1] - boundsa2j[0]; delta++)
      iintp2j[delta] = _li2.IntInterpolant(boundsa2j[0] + delta, ao2, bo2, _g2.GetJointGrid());

    // Get number of nodes of the joint grid along the first
    // direction.
    const int nxj1 = _g1.GetJointGrid().GetGrid().size();

    // Interpolate along the first direction of the joint grid
    std::vector<double> jg(nxj1, 0.);
    for (int beta = 0; beta < nxj1; beta++)
      for (int delta = 0; delta < boundsb2j[1] - boundsa2j[0]; delta++)
        jg[beta] += sgn2 * iintp2j[delta] * _dDJointGrid(beta, delta + boundsa2j[0]);

    // Now take care of the subgrids. ig2 has to correspond to the
    // subgrid index along the second direction such that the lower
    // bound "ao2" is on the denser grid possible.
    int ig2;
    for (ig2 = 0; ig2 < _g2.nGrids(); ig2++)
      if(_g2.GetSubGrid(ig2).xMin() > ao2)
        break;
    ig2--;

    // Get summation bounds on the igw-th subgrid along the second
    // direction.
    const std::array<int, 2> boundsa2s = _li2.SumBounds(ao2, _g2.GetSubGrid(ig2));
    const std::array<int, 2> boundsb2s = _li2.SumBounds(bo2, _g2.GetSubGrid(ig2));

    // Accumulate interpolating functions
    std::vector<double> iintp2s(boundsb2s[1] - boundsa2s[0]);
    for (int delta = 0; delta < boundsb2s[1] - boundsa2s[0]; delta++)
      iintp2s[delta] = _li2.IntInterpolant(boundsa2s[0] + delta, ao2, bo2, _g2.GetSubGrid(ig2));

    // Run over the subgrids along the first direction
    const int ng1 = _g1.nGrids();
    std::vector<std::vector<double>> sg(ng1);
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Get number of nodes of the ig1-th subgrid along the first
        // direction.
        const int nx1s = _g1.GetSubGrid(ig1).GetGrid().size();

        // Interpolate along the second direction of the ig1-th subgrid
        sg[ig1].resize(nx1s, 0.);
        for (int beta = 0; beta < nx1s; beta++)
          for (int delta = 0; delta < boundsb2s[1] - boundsa2s[0]; delta++)
            sg[ig1][beta] += sgn2 * iintp2s[delta] * _dDSubGrid[ig1][ig2](beta, delta + boundsa2j[0]);
      }
    return Distribution{_g1, sg, jg};
  }

  //_________________________________________________________________________
  DoubleDistribution DoubleDistribution::Derivative() const
  {
    return DoubleDistribution{this->_g1, this->_g2, [=] (double const& x1, double const& x2) -> double { return this->Derive(x1 - eps10, x2 - eps10); } };
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator = (DoubleDistribution const& d)
  {
    _dDSubGrid   = d._dDSubGrid;
    _dDJointGrid = d._dDJointGrid;

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator *= (double const& s)
  {
    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) *= s;

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
        for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
          for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
            _dDSubGrid[ig1][ig2](i, j) *= s;

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator *= (std::function<double(double const&, double const&)> const& f)
  {
    // Get joint grids
    const auto& jg1 = _g1.GetJointGrid().GetGrid();
    const auto& jg2 = _g2.GetJointGrid().GetGrid();

    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) *= f(jg1[i], jg2[j]);

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      {
        const auto& sg1 = _g1.GetSubGrid(ig1).GetGrid();
        for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
          {
            const auto& sg2 = _g2.GetSubGrid(ig2).GetGrid();
            for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
              for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
                _dDSubGrid[ig1][ig2](i, j) *= f(sg1[i], sg2[j]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator *= (std::function<double(double const&)> const& f)
  {
    // Get joint grids
    const auto& jg1 = _g1.GetJointGrid().GetGrid();
    const auto& jg2 = _g2.GetJointGrid().GetGrid();

    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) *= f(jg1[i]) * f(jg2[j]);

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      {
        const auto& sg1 = _g1.GetSubGrid(ig1).GetGrid();
        for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
          {
            const auto& sg2 = _g2.GetSubGrid(ig2).GetGrid();
            for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
              for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
                _dDSubGrid[ig1][ig2](i, j) *= f(sg1[i]) * f(sg2[j]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator /= (double const& s)
  {
    const double r = 1 / s;

    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) *= r;

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
        for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
          for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
            _dDSubGrid[ig1][ig2](i, j) *= r;

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator *= (DoubleDistribution const& d)
  {
    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) *= d._dDJointGrid(i, j);

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
        for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
          for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
            _dDSubGrid[ig1][ig2](i, j) *= _dDSubGrid[ig1][ig2](i, j);

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator += (DoubleDistribution const& d)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_g1 != &d._g1 || &this->_g2 != &d._g2)
      throw std::runtime_error(error("DoubleDistribution::operator +=", "DoubleDistribution grids do not match"));

    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) += d._dDJointGrid(i, j);

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
        for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
          for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
            _dDSubGrid[ig1][ig2](i, j) += _dDSubGrid[ig1][ig2](i, j);

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution& DoubleDistribution::operator -= (DoubleDistribution const& d)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_g1 != &d._g1 || &this->_g2 != &d._g2)
      throw std::runtime_error(error("DoubleDistribution::operator -=", "DoubleDistribution grids do not match"));

    // Joint grid
    for (size_t i = 0; i < _dDJointGrid.size(0); i++)
      for (size_t j = 0; j < _dDJointGrid.size(1); j++)
        _dDJointGrid(i, j) -= d._dDJointGrid(i, j);

    // Subgrids
    for (size_t ig1 = 0; ig1 < _dDSubGrid.size(); ig1++)
      for (size_t ig2 = 0; ig2 < _dDSubGrid[ig1].size(); ig2++)
        for (size_t i = 0; i < _dDSubGrid[ig1][ig2].size(0); i++)
          for (size_t j = 0; j < _dDSubGrid[ig1][ig2].size(1); j++)
            _dDSubGrid[ig1][ig2](i, j) -= _dDSubGrid[ig1][ig2](i, j);

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (double const& s, DoubleDistribution rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (DoubleDistribution lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (std::function<double(double const&, double const&)> const& f, DoubleDistribution rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (DoubleDistribution lhs, std::function<double(double const&, double const&)> const& f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (std::function<double(double const&)> const& f, DoubleDistribution rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (DoubleDistribution lhs, std::function<double(double const&)> const& f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  DoubleDistribution operator / (DoubleDistribution lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  DoubleDistribution operator + (DoubleDistribution lhs, DoubleDistribution const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  DoubleDistribution operator - (DoubleDistribution lhs, DoubleDistribution const& rhs)
  {
    return lhs -= rhs;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (DoubleDistribution lhs, DoubleDistribution const& rhs)
  {
    return lhs *= rhs;
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
