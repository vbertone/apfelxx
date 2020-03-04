//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/distribution.h"
#include "apfel/messages.h"

#include <stdexcept>
#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________
  Distribution::Distribution(Grid const& gr):
    LagrangeInterpolator{gr}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(Distribution                     const& obj,
                             std::vector<std::vector<double>> const& distsubgrid,
                             std::vector<double>              const& distjointgrid):
    LagrangeInterpolator{obj._grid}
  {
    _distributionSubGrid   = distsubgrid;
    _distributionJointGrid = distjointgrid;
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                             const& gr,
                             std::vector<std::vector<double>> const& distsubgrid,
                             std::vector<double>              const& distjointgrid):
    LagrangeInterpolator{gr}
  {
    _distributionSubGrid   = distsubgrid;
    _distributionJointGrid = distjointgrid;
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                 const& gr,
                             std::function<double(double const&)> const& InDistFunc):
    LagrangeInterpolator{gr}
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(InDistFunc(ix < 1 ? ix : 1));

    for (int ig = 0; ig < _grid.nGrids(); ig++)
      {
        std::vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          sg.push_back(InDistFunc(ix < 1 ? ix : 1));

        _distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                                const& gr,
                             std::function<double(double const&, double const&)> const& InDistFunc,
                             double                                              const& Q):
    LagrangeInterpolator{gr}
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(InDistFunc(ix < 1 ? ix : 1,Q));

    for (int ig = 0; ig < _grid.nGrids(); ig++)
      {
        std::vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          sg.push_back(InDistFunc(ix < 1 ? ix : 1,Q));

        _distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                             const& gr,
                             std::function<double(int const&, double const&)> const& InDistFunc,
                             int                                              const& ipdf):
    LagrangeInterpolator{gr}
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1));

    for (int ig = 0; ig < _grid.nGrids(); ig++)
      {
        std::vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          sg.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1));

        _distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                                            const& gr,
                             std::function<double(int const&, double const&, double const&)> const& InDistFunc,
                             int                                                             const& ipdf,
                             double                                                          const& Q):
    LagrangeInterpolator{gr}
  {
    for (auto const& ix: _grid.GetJointGrid().GetGrid())
      _distributionJointGrid.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1,Q));

    for (int ig = 0; ig < _grid.nGrids(); ig++)
      {
        std::vector<double> sg;
        for (auto const& ix: _grid.GetSubGrid(ig).GetGrid())
          sg.push_back(InDistFunc(ipdf,ix < 1 ? ix : 1,Q));

        _distributionSubGrid.push_back(sg);
      }
  }

  //_________________________________________________________________________
  void Distribution::PushJointGrid(double const& xi)
  {
    _distributionJointGrid.push_back(xi);
  }

  //_________________________________________________________________________
  void Distribution::PushSubGrid(double const& xi, bool const& next)
  {
    // If "next" is true start filling a new subgrid otherwise push
    // back in the last one.
    if (next)
      _distributionSubGrid.push_back(std::vector<double> {xi});
    else
      _distributionSubGrid.back().push_back(xi);
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
  Distribution& Distribution::operator *= (std::function<double(double const&)> const& f)
  {
    // Get joint grid
    const auto& jg = _grid.GetJointGrid().GetGrid();

    // sum objects in joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= f(jg[i]);

    // sum objects in subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      {
        // Get ig-th subgrid
        const auto& sg = _grid.GetSubGrid(ig).GetGrid();
        for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
          _distributionSubGrid[ig][i] *= f(sg[i]);
      }

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
      throw std::runtime_error(error("Distribution::operator+=", "Distribution grids does not match"));

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
      throw std::runtime_error(error("Distribution::operator+=", "Distribution grids does not match"));

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
  Distribution operator * (std::function<double(double const&)> const& f, Distribution rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  Distribution operator * (Distribution lhs, std::function<double(double const&)> const& f)
  {
    return lhs *= f;
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

  //_________________________________________________________________________
  Distribution operator * (Distribution lhs, Distribution const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  std::map<int,Distribution> DistributionMap(Grid                                                              const& g,
                                             std::function<std::map<int,double>(double const&, double const&)> const& InDistFunc,
                                             double                                                            const& Q,
                                             std::vector<int>                                                  const& skip)
  {
    // Joint grid and subgrid vectors.
    const std::vector<double>& jg = g.GetJointGrid().GetGrid();

    // Initialise output.
    std::map<int, Distribution> DistMap;
    const std::map<int, double> f = InDistFunc(jg[0], Q);
    for (auto it = f.begin(); it != f.end(); ++it)
      if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
        DistMap.insert({it->first, Distribution{g}});

    // Fill in joint grid.
    for (auto const& ix : jg)
      {
        const std::map<int, double> f = InDistFunc(ix < 1 ? ix : 1, Q);
        for (auto it = f.begin(); it != f.end(); ++it)
          if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
            DistMap.at(it->first).PushJointGrid(it->second);
      }

    // Fill in subgrids.
    for (int ig = 0; ig < g.nGrids(); ig++)
      {
        bool next = true;
        for (auto const& ix: g.GetSubGrid(ig).GetGrid())
          {
            const std::map<int, double> f = InDistFunc(ix < 1 ? ix : 1, Q);
            for (auto it = f.begin(); it != f.end(); ++it)
              if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
                DistMap.at(it->first).PushSubGrid(it->second, next);
            next = false;
          }
      }
    return DistMap;
  }

  //_________________________________________________________________________
  std::map<int,Distribution> DistributionMap(Grid                                               const& g,
                                             std::function<std::map<int,double>(double const&)> const& InDistFunc,
                                             std::vector<int>                                   const& skip)
  {
    // Joint grid and subgrid vectors.
    const std::vector<double>& jg = g.GetJointGrid().GetGrid();

    // Initialise output.
    std::map<int, Distribution> DistMap;
    const std::map<int, double> f = InDistFunc(jg[0]);
    for (auto it = f.begin(); it != f.end(); ++it)
      if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
        DistMap.insert({it->first, Distribution{g}});

    // Fill in joint grid.
    for (auto const& ix : jg)
      {
        const std::map<int, double> f = InDistFunc(ix < 1 ? ix : 1);
        for (auto it = f.begin(); it != f.end(); ++it)
          if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
            DistMap.at(it->first).PushJointGrid(it->second);
      }

    // Fill in subgrids.
    for (int ig = 0; ig < g.nGrids(); ig++)
      {
        bool next = true;
        for (auto const& ix: g.GetSubGrid(ig).GetGrid())
          {
            const std::map<int, double> f = InDistFunc(ix < 1 ? ix : 1);
            for (auto it = f.begin(); it != f.end(); ++it)
              if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
                DistMap.at(it->first).PushSubGrid(it->second, next);
            next = false;
          }
      }
    return DistMap;
  }
}
