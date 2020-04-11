//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/lagrangeinterpolator.h"
#include "apfel/constants.h"
#include "apfel/tools.h"

#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(Grid const& gr):
    Interpolator{gr}
  {
  }

  //_________________________________________________________________________________
  LagrangeInterpolator::LagrangeInterpolator(Grid                             const& gr,
                                             std::vector<std::vector<double>> const& distsubgrid,
                                             std::vector<double>              const& distjointgrid):
    Interpolator{gr, distsubgrid, distjointgrid}
  {
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::InterpolantLog(int const& beta, double const& lnx, SubGrid const& sg) const
  {
    // Get the logarithmic grid.
    const std::vector<double>& lxsg = sg.GetLogGrid();

    // Return immediately 1 if "x" coincides with "xg[beta]".
    if (std::abs(lnx - lxsg[beta]) < eps12)
      return 1;

    // Define the lower bound of the interpolation range.
    const int id    = sg.InterDegree();
    const int bound = std::max(beta - id, 0);

    // Return 0 if "x" is outside the range in which the interpolant
    // is different from zero.  Ideally this functions should never be
    // called if "beta" and "x" are such that "Interpolant" is
    // identically zero. Use "SumBounds" to know where "beta" should
    // run over given "x".
    if (lnx < lxsg[bound] || lnx >= lxsg[beta+1])
      return 0;

    // Find the the neighbors of "x" on the grid.
    int j;
    for (j = 0; j <= beta - bound; j++)
      if (lnx >= lxsg[beta-j])
        break;

    // Compute the interpolant.
    double w_int = 1;
    for (int delta = beta - j; delta <= beta - j + id; delta++)
      if (delta != beta)
        w_int *= ( lnx - lxsg[delta] ) / ( lxsg[beta] - lxsg[delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::Interpolant(int const& beta, double const& x, SubGrid const& sg) const
  {
    // Get the grid.
    const std::vector<double>& xg = sg.GetGrid();

    // Return immediately 1 if "x" coincides with "xg[beta]".
    if (std::abs(x - xg[beta]) < eps12)
      return 1;

    // Define the lower bound of the interpolation range.
    const int id    = sg.InterDegree();
    const int bound = std::max(beta - id, 0);

    // Return 0 if "x" is outside the range in which the interpolant
    // is different from zero.  Ideally this functions should never be
    // called if "beta" and "x" are such that "Interpolant" is
    // identically zero. Use "SumBounds" to know where "beta" should
    // run over given "x".
    if (x < xg[bound] || x >= xg[beta+1])
      return 0;

    // Find the the neighbors of "x" on the grid.
    int j;
    for (j = 0; j <= beta - bound; j++)
      if (x >= xg[beta-j])
        break;

    // Compute the interpolant.
    double w_int = 1;
    for (int delta = beta - j; delta <= beta - j + id; delta++)
      if (delta != beta)
        w_int *= ( x - xg[delta] ) / ( xg[beta] - xg[delta] );

    return w_int;
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::DerInterpolant(int const& beta, double const& x, SubGrid const& sg) const
  {
    // Get the grid.
    const std::vector<double>& xg = sg.GetGrid();

    // Define the lower bound of the interpolation range.
    const int id    = sg.InterDegree();
    const int bound = std::max(beta - id, 0);

    // Return 0 if "x" is outside the range in which the interpolant
    // is different from zero.  Ideally this functions should never be
    // called if "beta" and "x" are such that "Interpolant" is
    // identically zero. Use "SumBounds" to know where "beta" should
    // run over given "x".
    if (x < xg[bound] || x >= xg[beta+1])
      return 0;

    // Find the the neighbors of "x" on the grid.
    int j;
    for (j = 0; j <= beta - bound; j++)
      if (x >= xg[beta-j])
        break;

    // Compute the interpolant.
    double dw_int = 0;
    for (int gamma = beta - j; gamma <= beta - j + id; gamma++)
      {
        double w = 1;
        for (int delta = beta - j; delta <= beta - j + id; delta++)
          if (delta != beta && delta != gamma)
            w *= ( x - xg[delta] ) / ( xg[beta] - xg[delta] );
        if (gamma != beta)
          {
            w /= xg[beta] - xg[gamma];
            dw_int += w;
          }
      }
    return dw_int;
  }

  //_________________________________________________________________________________
  double LagrangeInterpolator::IntInterpolant(int const& beta, double const& a, double const& b, SubGrid const& sg) const
  {
    // Contrary to the other interpolants, here I use polynomials in
    // 'x' rather than in 'log(x)'. This seems to be faster and more
    // accurate. The algorithm in terms of polynomials in 'log(x)' is
    // appended below and commented out.

    // Get the grid.
    const std::vector<double>& xg = sg.GetGrid();

    // Interpolation degree
    const int k = sg.InterDegree();

    // Return 0 if "a" and "b" are outside the range in which the
    // interpolant is different from zero.
    if (a > xg[beta+1] || b < xg[std::max(beta-k, 0)])
      return 0;

    // Construct interpolant
    double iw_int = 0;
    for (int i = 0; i <= std::min(k, beta); i++)
      {
        if (xg[beta-i] > b || xg[beta-i+1] < a)
          continue;

        // Product of denominators
        double dp = 1;
        std::vector<double> r(k);

        int j = 0;
        for (int m = 0; m <= k; m++)
          if(m != i)
            {
              dp /= xg[beta] - xg[beta-i+m];
              r[j++] = xg[beta-i+m];
            }

        // Expansion coefficients
        const std::vector<double> p = ProductExpansion(r);

        // Integration bounds
        const double ab = std::max(a, xg[beta-i]);
        const double bb = std::min(b, xg[beta-i+1]);

        // Sum of the integrals
        double sum = 0;
        for (int n = 0; n <= k; n++)
          sum += pow(-1, n) * p[n] * ( pow(bb, k - n + 1) - pow(ab, k - n + 1) ) / ( k - n + 1 );

        iw_int += dp * sum;
      }
    return iw_int;
    /*
        // Let us first try using polynomials in 'log(x)'. This also works
        // but seems to me less accurate.
        // Get the grid.
        const std::vector<double>& lxg = sg.GetLogGrid();

        // Interpolation degree
        const int k = sg.InterDegree();

        // Transform bounds
        //const double lna = log(a);
        //const double lnb = log(b);

        // Construct interpolant
        double iw_int = 0;
        for (int i = 0; i <= std::min(k, beta); i++)
          {
            if (lxg[beta-i] > lnb || lxg[beta-i+1] < lna)
              continue;

            // Product of denominators
            double dp = 1;
            std::vector<double> r(k);
            int j = 0;
            for (int m = 0; m <= k; m++)
              if(m != i)
                {
                  dp /= lxg[beta] - lxg[beta-i+m];
                  r[j++] = lxg[beta-i+m];
                }

            // Expansion coefficients
            const std::vector<double> p = ProductExpansion(r);

            // Coefficients q
            std::vector<double> q(k+1, 0);
            for (int n = 0; n <= k; n++)
              for (int eta = n; eta <= k; eta++)
                q[n] += factorial(eta) * p[k-eta];

            // Integration bounds
            const double lab = std::max(lna, lxg[beta-i]);
            const double lbb = std::min(lnb, lxg[beta-i+1]);

            // Sum of the integrals
            double sum = 0;
            for (int n = 0; n <= k; n++)
              sum += q[n] * pow(-1, k + n) * ( exp(lbb) * pow(lbb, n) - exp(lab) * pow(lab, n) ) / factorial(n);

            iw_int += dp * sum;
          }
        return iw_int;
    */
  }

  //_________________________________________________________________________________
  std::array<int, 2> LagrangeInterpolator::SumBounds(double const& x, SubGrid const& sg) const
  {
    const std::vector<double>& xsg = sg.GetGrid();

    std::array<int,2> bounds = {{0, 0}};
    if (x < xsg[0] - eps12 || x > xsg[sg.nx()] + eps12)
      return bounds;

    const int low = std::lower_bound(xsg.begin(), xsg.end() - sg.InterDegree() - 1, x) - xsg.begin();
    bounds[0] = bounds[1] = low;

    if (std::abs(x - xsg[low]) <= eps12)
      bounds[1] += 1;
    else
      {
        bounds[0] -= 1;
        bounds[1] += sg.InterDegree();
      }

    return bounds;
  }
}
