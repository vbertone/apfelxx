//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/operatorgpd.h"
#include "apfel/integrator.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________
  OperatorGPD::OperatorGPD(Grid const& gr, Expression const& expr, double const& eps):
    Operator(gr)
  {
    // Interpolator object for the interpolating functions
    const LagrangeInterpolator li{_grid};

    // Loop over the subgrids
    _Operator.resize(1);

    // Get joint grid
    const SubGrid& jg = _grid.GetJointGrid();

    // Get vector with the grid nodes
    const std::vector<double>& xg = jg.GetGrid();

    // Get number of grid intervals in the definition range
    const int nx = jg.nx();

    // Get interpolation degree
    const int kappa = jg.InterDegree();

    // Initialise operator
    _Operator[0].resize(nx, nx, 0);
    // Loop over the index beta. In fact beta = 0 because the size
    // of the first dimension of "_Operator" is one.
    for (int beta = 0; beta < (int) _Operator[0].size(0); beta++)
      {
        // Set xg[beta] as external variable for the computation
        // of the operator (this is in general not needed but in
        // the case of GPD evolution).
        expr.SetExternalVariable(xg[beta]);

        // Loop over the index alpha
        for (int alpha = 0; alpha < (int) _Operator[0].size(1); alpha++)
          {
            // Weight of the subtraction term this is
            // \delta_{\alpha\beta}
            const double ws = (beta ==  alpha ? 1 : 0);

            // Given that the interpolation functions have
            // discontinuos derivative at the nodes, it turns out
            // that it is convenient to split the integrals into
            // kappa + 1 intervals on each of which the integrand is
            // smooth. Despite more integrals have to be computed,
            // the integration converges faster and is more
            // accurate.

            // Run over the grid intervals over which the
            // interpolating function is different from zero.
            for (int j = 0; j <= std::min(alpha, kappa); j++)
              {
                // Define "Integrator" object.
                const Integrator Ij{[&] (double const& y) -> double
                  {
                    const double wr = li.Interpolant(alpha, xg[beta] / y, jg);
                    return expr.Regular(y) * wr + expr.Singular(y) * ( wr - ws );
                  }};
                // Compute the integral
                _Operator[0](beta, alpha) += Ij.integrate(xg[beta] / xg[alpha - j + 1], xg[beta] / xg[alpha - j], eps);
              }
          }
        // Add the local part
        _Operator[0](beta, beta) += expr.Local(xg[beta] / xg[beta + 1]);
      }
  }

  //_________________________________________________________________________
  Distribution OperatorGPD::operator *= (Distribution const& d) const
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &d.GetGrid())
      throw std::runtime_error(error("OperatorGPD::operator *=", "Operator and Distribution grids do not match"));

    // Get number of subgrids
    const int ng = _grid.nGrids();

    // Get map of indices map from joint to subgrids
    const std::vector<std::vector<int>>& jsmap = _grid.JointToSubMap();

    // Get joint distribution
    const std::vector<double>& dj = d.GetDistributionJointGrid();

    // Initialise output vectors
    std::vector<double> j(d.GetDistributionJointGrid().size(), 0);
    std::vector<std::vector<double>> s(ng);

    // Get number of grid intervals in the definition range
    const int nx = _grid.GetJointGrid().nx();

    // Construct joint distribution first. The product between the
    // operator and the distribution is done exploiting the symmetry
    // of the operator.
    for (int beta = 0; beta < nx; beta++)
      for (int alpha = 0; alpha < nx; alpha++)
        j[beta] += _Operator[0](beta, alpha) * dj[alpha];

    // Compute the the distribution on the subgrids
    for (int ig = 0; ig < ng; ig++)
      {
        // Resize output on the "ig"-th subgrid
        s[ig].resize(d.GetDistributionSubGrid()[ig].size(), 0);
        for (int alpha = 0; alpha < _grid.GetSubGrid(ig).nx(); alpha++)
          s[ig][alpha] += j[jsmap[ig][alpha]];
      }

    // Return distribution object
    return Distribution{_grid, s, j};
  }

  //_________________________________________________________________________
  OperatorGPD& OperatorGPD::operator *= (OperatorGPD const& o)
  {
    // Fast method to check that we are using the same Grid
    if (&_grid != &o.GetGrid())
      throw std::runtime_error(error("OperatorGPD::operator*=", "Operators grid does not match"));

    const std::vector<matrix<double>> v = _Operator;
    // Set operator entries to zero
    _Operator[0].set(0);

    // The product between the operators is done exploiting the
    // symmetry of the operators.
    for (int alpha = 0; alpha < (int) _Operator[0].size(0); alpha++)
      for (int beta = 0; beta < (int) _Operator[0].size(1); beta++)
        for (int gamma = 0; gamma < (int) _Operator[0].size(1); gamma++)
          _Operator[0](alpha, beta) += v[0](alpha, gamma) * o._Operator[0](gamma, beta);

    return *this;
  }

  //_________________________________________________________________________
  Distribution operator * (OperatorGPD lhs, Distribution const& rhs)
  {
    return lhs *= rhs;
  }
}
