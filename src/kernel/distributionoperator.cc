//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/distributionoperator.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________
  DistributionOperator::DistributionOperator(Grid const& gr1, Grid const& gr2):
    _grid1(gr1),
    _grid2(gr2)
  {
    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Fill in matrix with zeros
    _dOperator.resize(ng1, std::vector<std::vector<matrix<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(nx1);
          for (int beta = 0; beta < nx1; beta++)
            _dOperator[ig1][ig2][beta].resize(1, nx2);
        }
  }

  //_________________________________________________________________________
  DistributionOperator::DistributionOperator(Distribution const& d1, Operator const& O2):
    _grid1(d1.GetGrid()),
    _grid2(O2.GetGrid())
  {
    // Stop the computation if the operators is of GPD type. Still
    // impossible to handle here.
    if (O2.IsGPD())
      throw std::runtime_error(error("DistributionOperator::DistributionOperator", "GPD operators cannot be handled yet."));

    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Get distribution and operator content
    const std::vector<std::vector<double>> d1g = d1.GetDistributionSubGrid();
    const std::vector<matrix<double>> O2g = O2.GetOperator();

    // Fill in matrix with zeros
    _dOperator.resize(ng1, std::vector<std::vector<matrix<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(nx1);
          for (int beta = 0; beta < (int) nx1; beta++)
            {
              _dOperator[ig1][ig2][beta].resize(1, nx2);
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(0); delta++)
                for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(1); gamma++)
                  _dOperator[ig1][ig2][beta](delta, gamma) = d1g[ig1][beta] * O2g[ig2](delta, gamma);
            }
        }
  }

  //_________________________________________________________________________
  DistributionOperator::DistributionOperator(DoubleObject<Distribution, Operator> const& DObj):
    _grid1(DObj.GetTerms()[0].object1.GetGrid()),
    _grid2(DObj.GetTerms()[0].object2.GetGrid())
  {
    // Stop the computation if the operator in the first term of the
    // double objects is of GPD type. Still impossible to handle.
    if (DObj.GetTerms()[0].object2.IsGPD())
      throw std::runtime_error(error("DistributionOperator::DistributionOperator", "GPD operators cannot be handled yet."));

    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Resize matrix and set to zero
    _dOperator.resize(ng1, std::vector<std::vector<matrix<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(nx1);
          for (int beta = 0; beta < nx1; beta++)
            _dOperator[ig1][ig2][beta].resize(1, nx2);
        }

    // Run over terms
    for (auto const& term : DObj.GetTerms())
      {
        // Get coefficient, distribution, and operator content of the
        // current term.
        const double coef = term.coefficient;
        const std::vector<std::vector<double>> d1g = term.object1.GetDistributionSubGrid();
        const std::vector<matrix<double>> O2g = term.object2.GetOperator();

        // Fill in matrix by multiplying distribution and operator of
        // the current term weighted by its coefficient.
        for (int ig1 = 0; ig1 < ng1; ig1++)
          for (int ig2 = 0; ig2 < ng2; ig2++)
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(0); delta++)
                for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(1); gamma++)
                  _dOperator[ig1][ig2][beta](delta, gamma) += coef * d1g[ig1][beta] * O2g[ig2](delta, gamma);
      }
  }

  //_________________________________________________________________________
  DoubleDistribution DistributionOperator::operator *= (Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid
    if (&_grid2 != &d.GetGrid())
      throw std::runtime_error(error("OperatorDistribution::operator *=", "DistributionOperator and Distribution grids do not match."));

    // Get numbers of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Get number of grid intervals in the definition range of the
    // joint grids.
    const int nx1 = _grid1.GetJointGrid().nx();
    const int nx2 = _grid2.GetJointGrid().nx();

    // Get maps of indices map from joint to subgrids
    const std::vector<std::vector<int>>& jsmap1 = _grid1.JointToSubMap();
    const std::vector<std::vector<int>>& jsmap2 = _grid2.JointToSubMap();

    // Get map of indices map from sub to joint to grids
    const std::vector<std::pair<int, int>>& sjmap1 = _grid1.SubToJointMap();
    const std::vector<std::pair<int, int>>& sjmap2 = _grid2.SubToJointMap();

    // Initialise output vectors
    matrix<double> j(_grid1.GetJointGrid().GetGrid().size(), _grid2.GetJointGrid().GetGrid().size());
    std::vector<std::vector<matrix<double>>> s(ng1, std::vector<matrix<double>>(ng2));

    // Get joint distribution
    const std::vector<double>& dj = d.GetDistributionJointGrid();

    // Construct joint distributions first. The product between the
    // operator and the distribution is done exploiting the symmetry
    // of the operator, which implies that it has one line
    // only for each pair of indices.
    for (int beta = 0; beta < nx1; beta++)
      {
        const std::pair<int, int> m1 = sjmap1[beta];
        for (int delta = 0; delta < nx2; delta++)
          {
            const std::pair<int, int> m2 = sjmap2[delta];
            for (int gamma = m2.second; gamma < _grid2.GetSubGrid(m2.first).nx(); gamma++)
              j(beta, delta) += _dOperator[m1.first][m2.first][m1.second](0, gamma - m2.second) * dj[jsmap2[m2.first][gamma]];
          }
      }

    // Compute the the distribution on the subgrids
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          // Resize output on the "ig1-ig2"-th subgrid
          s[ig1][ig2].resize(_grid1.GetSubGrid(ig1).GetGrid().size(), _grid2.GetSubGrid(ig2).GetGrid().size());
          for (int alpha = 0; alpha < _grid1.GetSubGrid(ig1).nx(); alpha++)
            for (int gamma = 0; gamma < _grid2.GetSubGrid(ig2).nx(); gamma++)
              s[ig1][ig2](alpha, gamma) += j(jsmap1[ig1][alpha], jsmap2[ig2][gamma]);
        }

    // Return double distribution object
    return DoubleDistribution{_grid1, _grid2, s, j};
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator *= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid2 != &o.GetGrid())
      throw std::runtime_error(error("DistributionOperator::operator *=", "Grids do not match."));

    // Get numbers of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Save current values
    const std::vector<std::vector<std::vector<matrix<double>>>> v = _dOperator;

    // Set current values to zero
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
          for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(0); delta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(1); gamma++)
              _dOperator[ig1][ig2][beta](delta, gamma) = 0;

    // Get operator values
    const std::vector<matrix<double>> ov = o.GetOperator();

    // Compute product
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        for (int beta = 0; beta < (int) v[ig1][ig2].size(); beta++)
          for (int delta = 0; delta < (int) v[ig1][ig2][beta].size(1); delta++)
            for (int sigma = 0; sigma <= delta; sigma++)
              _dOperator[ig1][ig2][beta](0, delta) += v[ig1][ig2][beta](0, sigma) * ov[ig2](0, delta - sigma);

    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator *= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
          for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(0); gamma++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(1); delta++)
              _dOperator[ig1][ig2][beta](gamma, delta) *= s;

    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator /= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
          for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(0); gamma++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(1); delta++)
              _dOperator[ig1][ig2][beta](gamma, delta) /= s;

    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator += (DistributionOperator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DistributionOperator::operator +=", "Grids do not match."));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
          for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(0); gamma++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(1); delta++)
              _dOperator[ig1][ig2][beta](gamma, delta) += o._dOperator[ig1][ig2][beta](gamma, delta);

    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator -= (DistributionOperator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DistributionOperator::operator -=", "Grids do not match."));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
          for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(0); gamma++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(1); delta++)
              _dOperator[ig1][ig2][beta](gamma, delta) -= o._dOperator[ig1][ig2][beta](gamma, delta);

    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator *= (std::function<double(double const&, double const&)> const& f)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      {
        const std::vector<double>& sg1 = _grid1.GetSubGrid(ig1).GetGrid();
        for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
          {
            const std::vector<double>& sg2 = _grid2.GetSubGrid(ig2).GetGrid();
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
              for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(0); gamma++)
                for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(1); delta++)
                  _dOperator[ig1][ig2][beta](gamma, delta) *= f(sg1[beta], sg2[delta]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator *= (std::function<double(double const&)> const& f)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      {
        const std::vector<double>& sg1 = _grid1.GetSubGrid(ig1).GetGrid();
        for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
          {
            const std::vector<double>& sg2 = _grid2.GetSubGrid(ig2).GetGrid();
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(); beta++)
              for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2][beta].size(0); gamma++)
                for (int delta = 0; delta < (int) _dOperator[ig1][ig2][beta].size(1); delta++)
                  _dOperator[ig1][ig2][beta](gamma, delta) *= f(sg1[beta]) * f(sg2[delta]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  DistributionOperator& DistributionOperator::operator = (DistributionOperator const& o)
  {
    if (this != &o)
      *this = o;

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (DistributionOperator lhs, Distribution const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (DistributionOperator lhs, Operator const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (double const& s, DistributionOperator rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (DistributionOperator lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  DistributionOperator operator / (DistributionOperator lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  DistributionOperator operator + (DistributionOperator lhs, DistributionOperator const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  DistributionOperator operator - (DistributionOperator lhs, DistributionOperator const& rhs)
  {
    return lhs -= rhs;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (std::function<double(double const&, double const&)> f, DistributionOperator rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (DistributionOperator lhs, std::function<double(double const&, double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (std::function<double(double const&)> f, DistributionOperator rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  DistributionOperator operator * (DistributionOperator lhs, std::function<double(double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, DistributionOperator const& dop)
  {
    const std::vector<std::vector<std::vector<matrix<double>>>> om = dop.GetDistributionOperator();
    os << "DistributionOperator: " << &dop << "\n";
    os << "DistributionOperator on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) om.size(); i++)
      for (int j = 0; j < (int) om[i].size(); j++)
        {
          os << "O[" << i << "][" << j << "]: [";
          for (int beta = 0; beta < (int) om[i][j].size(); beta++)
            for (int gamma = 0; gamma < (int) om[i][j][beta].size(0); gamma++)
              for (int delta = 0; delta < (int) om[i][j][beta].size(1); delta++)
                os << "{[" << beta << "](" << gamma << ", " << delta << ") : " << om[i][j][beta](gamma, delta) << "} ";
          os << "]\n";
        }
    os.copyfmt(default_format);
    return os;
  }
}
