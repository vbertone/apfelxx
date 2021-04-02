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
            const double ws = (beta == alpha ? 1 : 0);

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
                // Define "Integrator" object. IMPORTANT: the
                // particular form of the subtraction term only
                // applies to singular terms that behave as
                // 1/(1-y). If other singular functions are used, this
                // term has to be adjusted.
                const Integrator Ij{[&] (double const& y) -> double
                  {
                    const double wr = li.Interpolant(alpha, xg[beta] / y, jg);
                    return expr.Regular(y) * wr + expr.Singular(y) * ( wr - ws * ( 1 + (y > 1 ? ( 1 - y ) / y : 0) ) );
                  }};
                // Compute the integral
                _Operator[0](beta, alpha) += Ij.integrate(xg[beta] / xg[alpha - j + 1], xg[beta] / xg[alpha - j], eps);
              }
          }
        // Add the local parts: that from standard +-prescripted terms
        // ("Local") and that deriving from principal-valued integrals
        // ("LocalPV").
        _Operator[0](beta, beta) += expr.Local(xg[beta] / xg[beta + 1])
                                    + expr.LocalPV(xg[beta] / xg[beta + 1]) - (beta == 0 ? 0 : expr.LocalPV(xg[beta - std::min(beta, kappa)] / xg[beta]));
      }
  }
}
