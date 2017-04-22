//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/distribution.h"
#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________
  Distribution::Distribution(Grid const& gr):
    LagrangeInterpolator{gr}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(Distribution const& obj, vector<vector<double>> const& distsubgrid, vector<double> const& distjointgrid):
    LagrangeInterpolator{obj._grid}
  {
    _distributionSubGrid   = distsubgrid;
    _distributionJointGrid = distjointgrid;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator = (Distribution const& d)
  {
    _distributionSubGrid   = d.GetDistributionSubGrid();
    _distributionJointGrid = d.GetDistributionJointGrid();

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator *= (double const& s)
  {
    // sum objects in joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= s;

    // sum objects in subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] *= s;

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator /= (double const& s)
  {
    const double r = 1 / s;
    // sum objects in joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= r;

    // sum objects in subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] *= r;

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator *= (Distribution const& d)
  {
    // multiply objects in joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= d._distributionJointGrid[i];

    // sum objects in subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] *= d._distributionSubGrid[ig][i];

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator += (Distribution const& d)
  {
    // fast method to check that we are using the same Grid
    if (&this->_grid != &d._grid)
      throw runtime_exception("Distribution::operator+=", "Distribution grids does not match");

    // sum objects in joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] += d._distributionJointGrid[i];

    // sum objects in subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] += d._distributionSubGrid[ig][i];

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator -= (Distribution const& d)
  {
    // fast method to check that we are using the same Grid
    if (&this->_grid != &d._grid)
      throw runtime_exception("Distribution::operator+=", "Distribution grids does not match");

    // sum objects in joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] -= d._distributionJointGrid[i];

    // sum objects in subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] -= d._distributionSubGrid[ig][i];

    return *this;
  }

  //_________________________________________________________________________
  Distribution operator * (double const& s, Distribution rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  Distribution operator * (Distribution lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  Distribution operator / (Distribution lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  Distribution operator + (Distribution lhs, Distribution const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  Distribution operator - (Distribution lhs, Distribution const& rhs)
  {
    return lhs -= rhs;
  }
}
