//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/distribution.h"
#include "apfel/grid.h"
#include "apfel/subgrid.h"
#include "apfel/tools.h"

namespace apfel
{
  //_________________________________________________________________________
  Distribution::Distribution(Grid const& gr):
    LagrangeInterpolator{gr}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(const Distribution &obj, vector<vector<double>> const& distsubgrid, vector<double> const& distjointgrid):
    LagrangeInterpolator{obj._grid}
  {
    _distributionSubGrid   = distsubgrid;
    _distributionJointGrid = distjointgrid;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator=(Distribution const& rhs)
  {
    if(this != &rhs)
      {
	new(this) Distribution(rhs._grid);
	_distributionSubGrid   = rhs._distributionSubGrid;
	_distributionJointGrid = rhs._distributionJointGrid;
      }
    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator*=(double const& s)
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
  Distribution& Distribution::operator+=(Distribution const& d)
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
  Distribution& Distribution::operator-=(Distribution const& d)
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
  Distribution operator*(double const& s, Distribution rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  Distribution operator*(Distribution lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  Distribution operator+(Distribution lhs, Distribution const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  Distribution operator-(Distribution lhs, Distribution const& rhs)
  {
    return lhs -= rhs;
  }
}
