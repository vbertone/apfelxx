//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
    LagrangeInterpolator{obj._grid, distsubgrid, distjointgrid}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                             const& gr,
                             std::vector<std::vector<double>> const& distsubgrid,
                             std::vector<double>              const& distjointgrid):
    LagrangeInterpolator{gr, distsubgrid, distjointgrid}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                 const& gr,
                             std::function<double(double const&)> const& InDistFunc):
    LagrangeInterpolator{gr}
  {
    const std::vector<double>& jg = _grid.GetJointGrid().GetGrid();
    _distributionJointGrid.resize(jg.size());
    for (int ix = 0; ix < (int) jg.size(); ix++)
      _distributionJointGrid[ix] = InDistFunc(std::min(jg[ix], 1.));

    _distributionSubGrid.resize(_grid.nGrids());
    for (int ig = 0; ig < (int) _distributionSubGrid.size(); ig++)
      {
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        _distributionSubGrid[ig].resize(sg.size());
        for (int ix = 0; ix < (int) sg.size(); ix++)
          _distributionSubGrid[ig][ix] = InDistFunc(std::min(sg[ix], 1.));
      }
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                                const& gr,
                             std::function<double(double const&, double const&)> const& InDistFunc,
                             double                                              const& Q):
    LagrangeInterpolator{gr}
  {
    const std::vector<double>& jg = _grid.GetJointGrid().GetGrid();
    _distributionJointGrid.resize(jg.size());
    for (int ix = 0; ix < (int) jg.size(); ix++)
      _distributionJointGrid[ix] = InDistFunc(std::min(jg[ix], 1.), Q);

    _distributionSubGrid.resize(_grid.nGrids());
    for (int ig = 0; ig < (int) _distributionSubGrid.size(); ig++)
      {
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        _distributionSubGrid[ig].resize(sg.size());
        for (int ix = 0; ix < (int) sg.size(); ix++)
          _distributionSubGrid[ig][ix] = InDistFunc(std::min(sg[ix], 1.), Q);
      }
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                             const& gr,
                             std::function<double(int const&, double const&)> const& InDistFunc,
                             int                                              const& ipdf):
    LagrangeInterpolator{gr}
  {
    const std::vector<double>& jg = _grid.GetJointGrid().GetGrid();
    _distributionJointGrid.resize(jg.size());
    for (int ix = 0; ix < (int) jg.size(); ix++)
      _distributionJointGrid[ix] = InDistFunc(ipdf, std::min(jg[ix], 1.));

    _distributionSubGrid.resize(_grid.nGrids());
    for (int ig = 0; ig < (int) _distributionSubGrid.size(); ig++)
      {
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        _distributionSubGrid[ig].resize(sg.size());
        for (int ix = 0; ix < (int) sg.size(); ix++)
          _distributionSubGrid[ig][ix] = InDistFunc(ipdf, std::min(sg[ix], 1.));
      }
  }

  //_________________________________________________________________________
  Distribution::Distribution(Grid                                                            const& gr,
                             std::function<double(int const&, double const&, double const&)> const& InDistFunc,
                             int                                                             const& ipdf,
                             double                                                          const& Q):
    LagrangeInterpolator{gr}
  {
    const std::vector<double>& jg = _grid.GetJointGrid().GetGrid();
    _distributionJointGrid.resize(jg.size());
    for (int ix = 0; ix < (int) jg.size(); ix++)
      _distributionJointGrid[ix] = InDistFunc(ipdf, std::min(jg[ix], 1.), Q);

    _distributionSubGrid.resize(_grid.nGrids());
    for (int ig = 0; ig < (int) _distributionSubGrid.size(); ig++)
      {
        const std::vector<double>& sg = _grid.GetSubGrid(ig).GetGrid();
        _distributionSubGrid[ig].resize(sg.size());
        for (int ix = 0; ix < (int) sg.size(); ix++)
          _distributionSubGrid[ig][ix] = InDistFunc(ipdf, std::min(sg[ix], 1.), Q);
      }
  }

  //_________________________________________________________________________
  void Distribution::SetJointGrid(int const& ix, double const& x)
  {
    _distributionJointGrid[ix] = x;
  }

  //_________________________________________________________________________
  void Distribution::SetSubGrid(int const& ig, int const& ix, double const& x)
  {
    _distributionSubGrid[ig][ix] = x;
  }

  //_________________________________________________________________________
  Distribution Distribution::Derivative() const
  {
    return Distribution{this->_grid, [=] (double const& x) -> double { return this->Derive(x); } };
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
  std::map<int, Distribution> DistributionMap(Grid                                                               const& g,
                                              std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc,
                                              double                                                             const& Q,
                                              std::vector<int>                                                   const& skip)
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
    for (int ix = 0; ix < (int) jg.size(); ix++)
      {
        const std::map<int, double> f = InDistFunc(std::min(jg[ix], 1.), Q);
        for (auto const& it : f)
          if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
            DistMap.at(it.first).SetJointGrid(ix, it.second);
      }

    // Fill in subgrids.
    for (int ig = 0; ig < g.nGrids(); ig++)
      {
        const std::vector<double>& sg = g.GetSubGrid(ig).GetGrid();
        for (int ix = 0; ix < (int) sg.size(); ix++)
          {
            const std::map<int, double> f = InDistFunc(std::min(sg[ix], 1.), Q);
            for (auto const& it : f)
              if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
                DistMap.at(it.first).SetSubGrid(ig, ix, it.second);
          }
      }
    return DistMap;
  }

  //_________________________________________________________________________
  std::map<int, Distribution> DistributionMap(Grid                                                const& g,
                                              std::function<std::map<int, double>(double const&)> const& InDistFunc,
                                              std::vector<int>                                    const& skip)
  {
    // Joint grid and subgrid vectors.
    const std::vector<double>& jg = g.GetJointGrid().GetGrid();

    // Initialise output.
    std::map<int, Distribution> DistMap;
    const std::map<int, double> f = InDistFunc(jg[0]);
    for (auto const& it : f)
      if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
        DistMap.insert({it.first, Distribution{g}});

    // Fill in joint grid.
    for (int ix = 0; ix < (int) jg.size(); ix++)
      {
        const std::map<int, double> f = InDistFunc(std::min(jg[ix], 1.));
        for (auto const& it : f)
          if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
            DistMap.at(it.first).SetJointGrid(ix, it.second);
      }

    // Fill in subgrids.
    for (int ig = 0; ig < g.nGrids(); ig++)
      {
        const std::vector<double>& sg = g.GetSubGrid(ig).GetGrid();
        for (int ix = 0; ix < (int) sg.size(); ix++)
          {
            const std::map<int, double> f = InDistFunc(std::min(sg[ix], 1.));
            for (auto const& it : f)
              if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
                DistMap.at(it.first).SetSubGrid(ig, ix, it.second);
          }
      }
    return DistMap;
  }
}
