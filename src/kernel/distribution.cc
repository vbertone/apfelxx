//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/distribution.h"
#include "apfel/messages.h"
#include "apfel/constants.h"

#include <stdexcept>
#include <algorithm>
#include <numeric>

namespace apfel
{
  //_________________________________________________________________________
  Distribution::Distribution(Grid const& gr):
    LagrangeInterpolator{gr}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(Distribution const& obj):
    LagrangeInterpolator{obj._grid, obj._distributionSubGrid, obj._distributionJointGrid}
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
  void Distribution::SetJointGrid(std::vector<double> const& jg)
  {
    _distributionJointGrid = jg;
  }

  //_________________________________________________________________________
  void Distribution::SetSubGrid(int const& ig, std::vector<double> const& sg)
  {
    _distributionSubGrid[ig] = sg;
  }

  //_________________________________________________________________________
  void Distribution::SetSubGrids(std::vector<std::vector<double>> const& sgs)
  {
    _distributionSubGrid = sgs;
  }

  //_________________________________________________________________________
  Distribution Distribution::Derivative() const
  {
    return Distribution{this->_grid, [=] (double const& x) -> double { return this->Derive(x - eps10); } };
  }

  //_________________________________________________________________________
  Distribution Distribution::Transform(std::function<double(double const&)> const& TranformationFunc) const
  {
    // Transform joint grid
    std::vector<double> TdistributionJointGrid(_grid.GetJointGrid().GetGrid().size());
    std::transform(_distributionJointGrid.begin(), _distributionJointGrid.end(), TdistributionJointGrid.begin(), [=] (double const& f) -> double { return TranformationFunc(f); });

    // Transform subgrids
    std::vector<std::vector<double>> TdistributionSubGrid(_grid.nGrids());
    for (int ig = 0; ig < (int) TdistributionSubGrid.size(); ig++)
      {
        TdistributionSubGrid[ig].resize(_grid.GetSubGrid(ig).GetGrid().size());
        std::transform(_distributionSubGrid[ig].begin(), _distributionSubGrid[ig].end(), TdistributionSubGrid[ig].begin(), [=] (double const& f) -> double { return TranformationFunc(f); });
      }

    // Construct Distribution object and return
    return Distribution{_grid, TdistributionSubGrid, TdistributionJointGrid};
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator = (Distribution const& d)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &d._grid)
      throw std::runtime_error(error("Distribution::operator +=", "Distribution grids do not match"));

    _distributionSubGrid   = d.GetDistributionSubGrid();
    _distributionJointGrid = d.GetDistributionJointGrid();

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator *= (double const& s)
  {
    // Joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= s;

    // Subgrids
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

    // Joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= f(jg[i]);

    // Subgrids
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
    // Joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] /= s;

    // Subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] /= s;

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator *= (Distribution const& d)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &d._grid)
      throw std::runtime_error(error("Distribution::operator *=", "Distribution grids do not match"));

    // Joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] *= d._distributionJointGrid[i];

    // Subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] *= d._distributionSubGrid[ig][i];

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator += (Distribution const& d)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &d._grid)
      throw std::runtime_error(error("Distribution::operator +=", "Distribution grids do not match"));

    // Joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] += d._distributionJointGrid[i];

    // Subgrids
    for (size_t ig = 0; ig < _distributionSubGrid.size(); ig++)
      for (size_t i = 0; i < _distributionSubGrid[ig].size(); i++)
        _distributionSubGrid[ig][i] += d._distributionSubGrid[ig][i];

    return *this;
  }

  //_________________________________________________________________________
  Distribution& Distribution::operator -= (Distribution const& d)
  {
    // Fast method to check that we are using the same Grid
    if (&this->_grid != &d._grid)
      throw std::runtime_error(error("Distribution::operator -=", "Distribution grids do not match"));

    // Joint grid
    for (size_t i = 0; i < _distributionJointGrid.size(); i++)
      _distributionJointGrid[i] -= d._distributionJointGrid[i];

    // Subgrids
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
    // Joint grid and subgrid vectors
    const std::vector<double>& jg = g.GetJointGrid().GetGrid();

    // Initialise output
    std::map<int, Distribution> DistMap;
    const std::map<int, double> f = InDistFunc(jg[0], Q);
    for (auto it = f.begin(); it != f.end(); ++it)
      if (std::find(skip.begin(), skip.end(), it->first) == skip.end())
        DistMap.insert({it->first, Distribution{g}});

    // Fill in joint grid
    for (int ix = 0; ix < (int) jg.size(); ix++)
      {
        const std::map<int, double> f = InDistFunc(std::min(jg[ix], 1.), Q);
        for (auto const& it : f)
          if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
            DistMap.at(it.first).SetJointGrid(ix, it.second);
      }

    // Fill in subgrids
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
    // Joint grid and subgrid vectors
    const std::vector<double>& jg = g.GetJointGrid().GetGrid();

    // Initialise output
    std::map<int, Distribution> DistMap;
    const std::map<int, double> f = InDistFunc(jg[0]);
    for (auto const& it : f)
      if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
        DistMap.insert({it.first, Distribution{g}});

    // Fill in joint grid
    for (int ix = 0; ix < (int) jg.size(); ix++)
      {
        const std::map<int, double> f = InDistFunc(std::min(jg[ix], 1.));
        for (auto const& it : f)
          if (std::find(skip.begin(), skip.end(), it.first) == skip.end())
            DistMap.at(it.first).SetJointGrid(ix, it.second);
      }

    // Fill in subgrids
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

  //_________________________________________________________________________
  std::map<int, Distribution> DistributionMap(Grid                                              const& g,
                                              std::function<std::vector<double>(double const&)> const& InDistFunc,
                                              int                                               const& NOutputs)
  {
    // Joint grid and subgrid vectors
    const std::vector<double>& jg = g.GetJointGrid().GetGrid();

    // Initialise output. If the number of outputs is provided use
    // that, otherwise call the function at the first grid point.
    const int n = ( NOutputs == 0 ? InDistFunc(jg[0]).size() : NOutputs );
    std::map<int, Distribution> DistMap;
    for (int i = 0; i < n; i++)
      DistMap.insert({i, Distribution{g}});

    // Fill in joint grid
    for (int ix = 0; ix < (int) jg.size(); ix++)
      {
        const std::vector<double> f = InDistFunc(std::min(jg[ix], 1.));
        int i = 0;
        for (double const& v : f)
          DistMap.at(i++).SetJointGrid(ix, v);
      }

    // Fill in subgrids. Since the subgrids are locked the joint grid
    // already contains all the nodes therefore there is no need to
    // call the input function again.
    const std::vector<std::vector<int>>& m = g.JointToSubMap();
    for (int o = 0; o < n; o++)
      {
        const std::vector<double>& jv = DistMap.at(o).GetDistributionJointGrid();
        for (int ig = 0; ig < (int) m.size(); ig++)
          for (int ix = 0; ix < (int) m[ig].size(); ix++)
            DistMap.at(o).SetSubGrid(ig, ix, jv[m[ig][ix]]);
      }
    return DistMap;
  }

  //_________________________________________________________________________
  double Sum(Distribution const& InDist)
  {
    const std::vector<double>& jv = InDist.GetDistributionJointGrid();
    return std::accumulate(jv.begin(), jv.end(), 0.);
  }

  //_________________________________________________________________________
  double InnerProduct(Distribution const& d1, Distribution const& d2, double const& offset)
  {
    return Sum(d1 * d2) + offset;
  }
}
