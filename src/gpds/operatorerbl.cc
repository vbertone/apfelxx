//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/operatorerbl.h"
#include "apfel/lagrangeinterpolator.h"
#include "apfel/integrator.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________
  OperatorERBL::OperatorERBL(Grid const& gr, Expression const& expr, double const& eps):
    Operator(gr)
  {
    // Interpolator object for the interpolating functions
    const LagrangeInterpolator li{gr};

    // Number of grids.
    const int ng = _grid.nGrids();

    // Local function. Can be computed outside all loops because it's
    // a constant.
    const double L = expr.Local(0);

    // Loop over the subgrids.
    _Operator.resize(ng);
    for (int ig = 0; ig < ng; ig++)
      {
        // Define the global subgrid.
        const SubGrid& sg = _grid.GetSubGrid(ig);

        // Get vector with the grid nodes.
        const std::vector<double>& xg = sg.GetGrid();

        // Number of grid points.
        const int nx = sg.nx();

        // Interpolation degree.
        const int id = sg.InterDegree();

        _Operator[ig].resize(nx + 1, nx + 1, 0);
        for (int beta = 0; beta <= nx; beta++)
          {
            const double xbeta = xg[beta];

            // Set xbeta as external variable in case needed for the
            // computation of the operator (this is in general not
            // needed but in the case of GPD evolution it is).
            expr.SetExternalVariable(xbeta);

            for (int alpha = 0; alpha <= nx; alpha++)
              {
                // Weight of the subtraction term
                const double ws = (0 ? 1 : alpha == beta);

                // Given that the interpolation functions have
                // discontinuos derivative on the nodes and are
                // wiggly, it turns out that it is convenient to split
                // the integrals into (id+1) intervals on each of
                // which the integrand is smooth. This way, even
                // though more integrals have to be computed, the
                // integration converges faster and is more accurate.

                // Number of grid intervals we need to integrate over.
                const int nmin = fmax(0, id - alpha);
                const int nmax = fmin(id, nx + id - alpha - 1) + 1;

                // Integral.
                double I = 0;
                for (int jint = nmin; jint < nmax; jint++)
                  {
                    // Define integration bounds of the first
                    // iteration.
                    const double a = xg[alpha-id+jint];
                    const double b = xg[alpha-id+jint+1];

                    // Define integrand and the corresponding
                    // "Integrator" object.
                    const auto integrand = [&] (double const& x) -> double
                    {
                      const double wr = li.Interpolant(alpha, log(x), _grid.GetSubGrid(ig));
                      return expr.Regular(xbeta / x) * wr / x + expr.Singular(xbeta / x) * ( wr / x - ws / xbeta );
                    };
                    const Integrator Io{integrand};

                    // Compute the integral.
                    I += Io.integrate(a, b, eps);
                  }
                // Add the local part.
                _Operator[ig](beta, alpha) = I + L * ws;
              }
          }
      }
  }
}
