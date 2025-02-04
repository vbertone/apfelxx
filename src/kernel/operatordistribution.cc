//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/operatordistribution.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________
  OperatorDistribution::OperatorDistribution(Grid const& gr1, Grid const& gr2):
    _grid1(gr1),
    _grid2(gr2)
  {
    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Fill in matrix with zeros
    _dOperator.resize(ng1, std::vector<matrix<std::vector<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(1, nx1, {});
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
            for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
              _dOperator[ig1][ig2](beta, alpha).resize(nx2);
        }
  }

  //_________________________________________________________________________
  OperatorDistribution::OperatorDistribution(Grid const& gr1, Grid const& gr2, std::vector<std::vector<matrix<std::vector<double>>>> const& OD):
    _grid1(gr1),
    _grid2(gr2),
    _dOperator(OD)
  {
    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Make sure that grid sizes match
    if ((int) OD.size() != ng1 || (int) OD[0].size() != ng2)
      throw std::runtime_error(error("OperatorDistribution::OperatorDistribution", "Numbers of subgrids of input grids and DistributionOperator object do not match"));

    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          if ((int) OD[ig1][ig2].size(1) != _grid1.GetSubGrid(ig1).nx())
            throw std::runtime_error(error("OperatorDistribution::OperatorDistribution", "Numbers of subgrid nodes of first grid do not match"));

          if ((int) OD[ig1][ig2](0, 0).size() != _grid2.GetSubGrid(ig2).nx())
            throw std::runtime_error(error("OperatorDistribution::OperatorDistribution", "Numbers of subgrid nodes of second grid do not match"));
        }
  }

  //_________________________________________________________________________
  OperatorDistribution::OperatorDistribution(Operator const& O1, Distribution const& d2):
    _grid1(O1.GetGrid()),
    _grid2(d2.GetGrid())
  {
    // Stop the computation if the operators is of GPD type. Still
    // impossible to handle here.
    if (O1.IsGPD())
      throw std::runtime_error(error("OperatorDistribution::OperatorDistribution", "GPD operators cannot be handled yet."));

    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Get distribution and operator content
    const std::vector<matrix<double>> O1g = O1.GetOperator();
    const std::vector<std::vector<double>> d2g = d2.GetDistributionSubGrid();

    // Fill in matrix with zeros
    _dOperator.resize(ng1, std::vector<matrix<std::vector<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(1, nx1, {});
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
            for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
              {
                _dOperator[ig1][ig2](beta, alpha).resize(nx2);
                for (int delta = 0; delta < (int) nx1; delta++)
                  _dOperator[ig1][ig2](beta, alpha)[delta] = O1g[ig1](beta, alpha) * d2g[ig2][delta];
              }
        }
  }

  //_________________________________________________________________________
  OperatorDistribution::OperatorDistribution(DoubleObject<Operator, Distribution> const& DObj):
    _grid1(DObj.GetTerms()[0].object1.GetGrid()),
    _grid2(DObj.GetTerms()[0].object2.GetGrid())
  {
    // Stop the computation if the operator in the first term of the
    // double objects is of GPD type. Still impossible to handle.
    if (DObj.GetTerms()[0].object1.IsGPD())
      throw std::runtime_error(error("OperatorDistribution::OperatorDistribution", "GPD operators cannot be handled yet."));

    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Resize matrix and set to zero
    _dOperator.resize(ng1, std::vector<matrix<std::vector<double>>>(ng2));
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        {
          const int nx1 = _grid1.GetSubGrid(ig1).nx();
          const int nx2 = _grid2.GetSubGrid(ig2).nx();
          _dOperator[ig1][ig2].resize(1, nx1, {});
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
            for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
              _dOperator[ig1][ig2](beta, alpha).resize(nx2);
        }

    // Run over terms
    for (auto const& term : DObj.GetTerms())
      {
        // Get coefficient, operator, and distribution content of the
        // current term.
        const double coef = term.coefficient;
        const std::vector<matrix<double>> O1g = term.object1.GetOperator();
        const std::vector<std::vector<double>> d2g = term.object2.GetDistributionSubGrid();

        // Fill in matrix by multiplying distribution and operator of
        // the current term weighted by its coefficient.
        for (int ig1 = 0; ig1 < ng1; ig1++)
          for (int ig2 = 0; ig2 < ng2; ig2++)
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
                  _dOperator[ig1][ig2](beta, alpha)[delta] += coef * O1g[ig1](beta, alpha) * d2g[ig2][delta];
      }
  }

  //_________________________________________________________________________
  DoubleDistribution OperatorDistribution::operator *= (Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &d.GetGrid())
      throw std::runtime_error(error("OperatorDistribution::operator *=", "OperatorDistribution and Distribution grids do not match."));

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
            for (int alpha = m1.second; alpha < _grid1.GetSubGrid(m1.first).nx(); alpha++)
              j(beta, delta) += _dOperator[m1.first][m2.first](0, alpha - m2.second)[m2.second] * dj[jsmap1[m1.first][alpha]];
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
  OperatorDistribution& OperatorDistribution::operator *= (Operator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetGrid())
      throw std::runtime_error(error("OperatorDistribution::operator *=", "Grids do not match."));

    // Get numbers of subgrids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Save current values
    const std::vector<std::vector<matrix<std::vector<double>>>> v = _dOperator;

    // Set current values to zero
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
          for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
              _dOperator[ig1][ig2](beta, alpha)[delta] = 0;

    // Get operator values
    const std::vector<matrix<double>> ov = o.GetOperator();

    // Compute product
    for (int ig1 = 0; ig1 < ng1; ig1++)
      for (int ig2 = 0; ig2 < ng2; ig2++)
        for (int alpha = 0; alpha < (int) v[ig1][ig2].size(1); alpha++)
          for (int delta = 0; delta < (int) v[ig1][ig2](0, alpha).size(); delta++)
            for (int rho = 0; rho <= alpha; rho++)
              _dOperator[ig1][ig2](0, alpha)[delta] += v[ig1][ig2](0, rho)[delta] * ov[ig1](0, alpha - rho);

    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator *= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
          for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
              _dOperator[ig1][ig2](beta, alpha)[delta] *= s;

    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator /= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
          for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
              _dOperator[ig1][ig2](beta, alpha)[delta] /= s;

    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator += (OperatorDistribution const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("OperatorDistribution::operator +=", "Grids do not match."));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
          for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
              _dOperator[ig1][ig2](beta, alpha)[delta] += o._dOperator[ig1][ig2](beta, alpha)[delta];

    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator -= (OperatorDistribution const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("OperatorDistribution::operator -=", "Grids do not match."));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
          for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
            for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
              _dOperator[ig1][ig2](beta, alpha)[delta] -= o._dOperator[ig1][ig2](beta, alpha)[delta];

    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator *= (std::function<double(double const&, double const&)> const& f)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      {
        const std::vector<double>& sg1 = _grid1.GetSubGrid(ig1).GetGrid();
        for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
          {
            const std::vector<double>& sg2 = _grid2.GetSubGrid(ig2).GetGrid();
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
                  _dOperator[ig1][ig2](beta, alpha)[delta] *= f(sg1[alpha], sg2[delta]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator *= (std::function<double(double const&)> const& f)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      {
        const std::vector<double>& sg1 = _grid1.GetSubGrid(ig1).GetGrid();
        for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
          {
            const std::vector<double>& sg2 = _grid2.GetSubGrid(ig2).GetGrid();
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                for (int delta = 0; delta < (int) _dOperator[ig1][ig2](beta, alpha).size(); delta++)
                  _dOperator[ig1][ig2](beta, alpha)[delta] *= f(sg1[alpha]) * f(sg2[delta]);
          }
      }
    return *this;
  }

  //_________________________________________________________________________
  OperatorDistribution& OperatorDistribution::operator = (OperatorDistribution const& o)
  {
    if (this != &o)
      *this = o;

    return *this;
  }

  //_________________________________________________________________________
  DoubleDistribution operator * (OperatorDistribution const& lhs, Distribution const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (OperatorDistribution lhs, Operator const& rhs)
  {
    return lhs *= rhs;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (double const& s, OperatorDistribution rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (OperatorDistribution lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________
  OperatorDistribution operator / (OperatorDistribution lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________
  OperatorDistribution operator + (OperatorDistribution lhs, OperatorDistribution const& rhs)
  {
    return lhs += rhs;
  }

  //_________________________________________________________________________
  OperatorDistribution operator - (OperatorDistribution lhs, OperatorDistribution const& rhs)
  {
    return lhs -= rhs;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (std::function<double(double const&, double const&)> f, OperatorDistribution rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (OperatorDistribution lhs, std::function<double(double const&, double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (std::function<double(double const&)> f, OperatorDistribution rhs)
  {
    return rhs *= f;
  }

  //_________________________________________________________________________
  OperatorDistribution operator * (OperatorDistribution lhs, std::function<double(double const&)> f)
  {
    return lhs *= f;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, OperatorDistribution const& dop)
  {
    const std::vector<std::vector<matrix<std::vector<double>>>> om = dop.GetOperatorDistribution();
    os << "OperatorDistribution: " << &dop << "\n";
    os << "OperatorDistribution on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) om.size(); i++)
      for (int j = 0; j < (int) om[i].size(); j++)
        {
          os << "O[" << i << "][" << j << "]: [";
          for (int alpha = 0; alpha < (int) om[i][j].size(0); alpha++)
            for (int beta = 0; beta < (int) om[i][j].size(1); beta++)
              for (int delta = 0; delta < (int) om[i][j](alpha, beta).size(); delta++)
                os << "{(" << alpha << ", " << beta << ")[" << delta << "] : " << om[i][j](alpha, beta)[delta] << "} ";
          os << "]\n";
        }
    os.copyfmt(default_format);
    return os;
  }
}
