//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubleoperator.h"
#include "apfel/operator.h"
#include "apfel/integrator.h"
#include "apfel/messages.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________
  DoubleOperator::DoubleOperator(Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr, double const& eps):
    _grid1(gr1),
    _grid2(gr2),
    _dexpr(dexpr),
    _eps(eps)
  {
    // Interpolator objects for the interpolating functions
    const LagrangeInterpolator li1{_grid1};
    const LagrangeInterpolator li2{_grid2};

    // Number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Resize operator w.r.t. the first grid
    _dOperator.resize(ng1);

    // Loop over the subgrids of the first grid
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        // Get current subgrid parameters on the first grid
        const SubGrid& sg1 = _grid1.GetSubGrid(ig1);
        const std::vector<double>& xg1 = sg1.GetGrid();
        const int nx1 = sg1.nx();
        const int id1 = sg1.InterDegree();
        const double dx1 = sg1.Step();

        // Resize operator w.r.t. the second grid
        _dOperator[ig1].resize(ng2);

        // Loop over the subgrids of the second grid
        for (int ig2 = 0; ig2 < ng2; ig2++)
          {
            // Get current subgrid parameters on the second grid
            const SubGrid& sg2 = _grid2.GetSubGrid(ig2);
            const std::vector<double>& xg2 = sg2.GetGrid();
            const int nx2 = sg2.nx();
            const int id2 = sg2.InterDegree();
            const double dx2 = sg2.Step();

            // Evaluate Local-Local contribution
            const double cLL = _dexpr.LocalLocal(1 / exp(dx1), 1 / exp(dx2));

            // Evaluate Local-Singular and Local-Regular on the second grid
            matrix<double> cLS_LR{1, (size_t) nx2};
            for (int delta = 0; delta < (int) cLS_LR.size(0); delta++)
              for (int gamma = delta; gamma < (int) cLS_LR.size(1); gamma++)
                for (int k = 0; k < std::min(id2, gamma - delta) + 1; k++)
                  {
                    const double ws2 = (delta == gamma ? 1 : 0);
                    const Integrator Ik{[&] (double const& y2) -> double
                      {
                        const double wr2 = li2.InterpolantLog(gamma, log(xg2[delta] / y2), sg2);
                        return _dexpr.LocalSingular(1 / exp(dx1), y2) * ( wr2 - ws2 ) + _dexpr.LocalRegular(1 / exp(dx1), y2) * wr2;
                      }};
                    cLS_LR(delta, gamma) += Ik.integrate(exp((delta - gamma + k - 1) * dx2), exp((delta - gamma + k) * dx2), _eps);
                  }

            // Evaluate Singular-Local and Regular-Local on the first grid
            matrix<double> cSL_RL{1, (size_t) nx1};
            for (int beta = 0; beta < (int) cSL_RL.size(0); beta++)
              for (int alpha = beta; alpha < (int) cSL_RL.size(1); alpha++)
                for (int j = 0; j < std::min(id1, alpha - beta) + 1; j++)
                  {
                    const double ws1 = (beta == alpha ? 1 : 0);
                    const Integrator Ij{[&] (double const& y1) -> double
                      {
                        const double wr1 = li1.InterpolantLog(alpha, log(xg1[beta] / y1), sg1);
                        return ( wr1 - ws1 ) * _dexpr.SingularLocal(y1, 1 / exp(dx2)) + wr1 * _dexpr.RegularLocal(y1, 1 / exp(dx2));
                      }};
                    cSL_RL(beta, alpha) += Ij.integrate(exp((beta - alpha + j - 1) * dx1), exp((beta - alpha + j) * dx1), _eps);
                  }

            // Finally move to compute Singular-Singular,
            // Singular-Regular, Regular-Singular, and Regular-Regular
            // that require integrating over both y1 and y2. Include
            // Local-Local, Local-Singular, Local-Regular,
            // Singular-Local, and Regular-Local along the way.

            // Resize operator w.r.t. the first variable
            _dOperator[ig1][ig2].resize(1, nx1);

            // Loop over the index beta. In fact beta = 0 because the size
            // of the first dimension of "_Operator" is one.
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              {
                // Loop over the index alpha
                for (int alpha = beta; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                  {
                    // Resize operator w.r.t. the second variable
                    _dOperator[ig1][ig2](beta, alpha).resize(1, nx2);

                    // Weight of the subtraction term w.r.t the first grid
                    const double ws1 = (beta == alpha ? 1 : 0);

                    // Loop over the index delta. In fact delta = 0 because the size
                    // of the first dimension of "_Operator" is one.
                    for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(0); delta++)
                      {
                        // Loop over the index gamma
                        for (int gamma = delta; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(1); gamma++)
                          {
                            // Weight of the subtraction term w.r.t the first grid
                            const double ws2 = (delta == gamma ? 1 : 0);

                            // Initialise integral
                            double cSS_SR_RS_RR = 0;

                            // Run over the grid intervals
                            for (int j = 0; j < std::min(id1, alpha - beta) + 1; j++)
                              {
                                for (int k = 0; k < std::min(id2, gamma - delta) + 1; k++)
                                  {
                                    const Integrator Iy1{[&] (double const& y1) -> double
                                      {
                                        const double wr1 = li1.InterpolantLog(alpha, log(xg1[beta] / y1), sg1);
                                        const Integrator Iy1y2{
                                          [&] (double const& y2) -> double
                                          {
                                            const double wr2 = li2.InterpolantLog(gamma, log(xg2[delta] / y2), sg2);
                                            return ( wr1 - ws1 ) * _dexpr.SingularSingular(y1, y2) * ( wr2 - ws2 )
                                                                         + ( wr1 - ws1 ) * _dexpr.SingularRegular(y1, y2) * wr2
                                                                         + wr1 * _dexpr.RegularSingular(y1, y2) * ( wr2 - ws2 )
                                                                         + wr1 * _dexpr.RegularRegular(y1, y2) * wr2;
                                          }};
                                        return Iy1y2.integrate(exp((delta - gamma + k - 1) * dx2), exp((delta - gamma + k) * dx2), _eps);
                                      }};
                                    cSS_SR_RS_RR += Iy1.integrate(exp((beta - alpha + j - 1) * dx1), exp((beta - alpha + j) * dx1), _eps);
                                  }
                              }
                            _dOperator[ig1][ig2](beta, alpha)(delta, gamma) = ws1 * cLL * ws2 + ws1 * cLS_LR(delta, gamma) + cSL_RL(beta, alpha) * ws2 + cSS_SR_RS_RR;
                          }
                      }
                  }
              }
          }
      }
  }

  //_________________________________________________________________________
  DoubleOperator::DoubleOperator(Operator const& O1, Operator const& O2):
    _grid1(O1.GetGrid()),
    _grid2(O2.GetGrid()),
    _dexpr(DoubleExpression{O1.GetExpression(), O2.GetExpression()}),
    _eps(std::max(O1.GetIntegrationAccuracy(), O2.GetIntegrationAccuracy()))
  {
    // Get single operators
    const std::vector<matrix<double>> o1 = O1.GetOperator();
    const std::vector<matrix<double>> o2 = O2.GetOperator();

    // Get number of grids
    const int ng1 = _grid1.nGrids();
    const int ng2 = _grid2.nGrids();

    // Fill in double operator
    _dOperator.resize(ng1);
    for (int ig1 = 0; ig1 < ng1; ig1++)
      {
        const int nx1 = _grid1.GetSubGrid(ig1).nx();
        _dOperator[ig1].resize(ng2);
        for (int ig2 = 0; ig2 < ng2; ig2++)
          {
            const int nx2 = _grid2.GetSubGrid(ig2).nx();
            _dOperator[ig1][ig2].resize(1, nx1);
            for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(0); beta++)
              for (int alpha = beta; alpha < (int) _dOperator[ig1][ig2].size(1); alpha++)
                {
                  _dOperator[ig1][ig2](beta, alpha).resize(1, nx2);
                  for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(0); delta++)
                    for (int gamma = delta; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(1); gamma++)
                      _dOperator[ig1][ig2](beta, alpha)(delta, gamma) = o1[ig1](beta, alpha) * o2[ig2](delta, gamma);
                }
          }
      }
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator *= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) *= s;

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator /= (double const& s)
  {
    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) /= s;

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator += (DoubleOperator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DoubleOperator::operator +=", "Grids do not match"));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) += o._dOperator[ig1][ig2](alpha, beta)(gamma, delta);

    return *this;
  }

  //_________________________________________________________________________
  DoubleOperator& DoubleOperator::operator -= (DoubleOperator const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid1 != &o.GetFirstGrid() || &_grid2 != &o.GetSecondGrid())
      throw std::runtime_error(error("DoubleOperator::operator +=", "Grids do not match"));

    for (int ig1 = 0; ig1 < (int) _dOperator.size(); ig1++)
      for (int ig2 = 0; ig2 < (int) _dOperator[ig1].size(); ig2++)
        for (int alpha = 0; alpha < (int) _dOperator[ig1][ig2].size(0); alpha++)
          for (int beta = 0; beta < (int) _dOperator[ig1][ig2].size(1); beta++)
            for (int gamma = 0; gamma < (int) _dOperator[ig1][ig2](alpha, beta).size(0); gamma++)
              for (int delta = 0; delta < (int) _dOperator[ig1][ig2](alpha, beta).size(1); delta++)
                _dOperator[ig1][ig2](alpha, beta)(gamma, delta) -= o._dOperator[ig1][ig2](alpha, beta)(gamma, delta);

    return *this;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, DoubleOperator const& dop)
  {
    const std::vector<std::vector<matrix<matrix<double>>>> om = dop.GetDoubleOperator();
    os << "Operator: " << &dop << "\n";
    os << "Operator on the SubGrids:" << "\n";
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(2);
    for (int i = 0; i < (int) om.size(); i++)
      for (int j = 0; j < (int) om[i].size(); j++)
        {
          os << "O[" << i << "][" << j << "]: [";
          for (int alpha = 0; alpha < (int) om[i][j].size(0); alpha++)
            for (int beta = 0; beta < (int) om[i][j].size(1); beta++)
              for (int gamma = 0; gamma < (int) om[i][j](alpha, beta).size(0); gamma++)
                for (int delta = 0; delta < (int) om[i][j](alpha, beta).size(1); delta++)
                  os << "{(" << alpha << ", " << beta << ")(" << gamma << ", " << delta << ") : " << om[i][j](alpha, beta)(gamma, delta) << "} ";
          os << "]\n";
        }
    os.copyfmt(default_format);
    return os;
  }
}
