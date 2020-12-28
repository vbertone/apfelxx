//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/integrator.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________
  Integrator::Integrator(std::function<double(double const&)> const& func, IntegrationMethod const& method):
    _func(func),
    _method(method)
  {
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, double const& eps) const
  {
    switch (_method)
      {
      case GAUSS_LEGENDRE:
        return integrateGL(a, b, eps);
      case GAUSS_KRONROD:
        return integrateGK(a, b, eps);
      default:
        throw std::runtime_error(error("Integrator::integrate", "Unknown integration method"));
      }
  }

  //_________________________________________________________________________
  double Integrator::integrateGL(double const& a, double const& b, double const& eps) const
  {
    const double delta = eps15 * std::abs(a - b);
    double dgauss = 0;
    double aa = a;

    while (true)
      {
        double y = b - aa;
        if (std::abs(y) < delta)
          break;

        double s1;
        double s0;
        do
          {
            const double bb = aa + y;

            // Central point and half interval
            const double c1 = ( bb + aa ) / 2;
            const double c2 = ( bb - aa ) / 2;

            // Integral with smaller number of points
            s0 = 0;
            for (int i = 0; i < (int) gl_x[0].size(); i++)
              {
                const double u = gl_x[0][i] * c2;
                s0 += gl_w[0][i] * ( _func(c1 + u) + _func(c1 - u) );
              }
            s0 *= c2;

            // Integral with larger number of points
            s1 = 0;
            for (int i = 0; i < (int) gl_x[1].size(); i++)
              {
                const double u = gl_x[1][i] * c2;
                s1 += gl_w[1][i] * ( _func(c1 + u) + _func(c1 - u) );
              }
            s1 *= c2;

            // Split the interval into two
            y /= 2;

            // Stop if the interval is smaller than the minimum step
            if (std::abs(y) < delta)
              throw std::runtime_error(error("Integrator::integrateGL", "Too high accuracy required."));

            // Check accuracy and continue if necessary
          }
        while (std::abs(s1 - s0) / ( 1 + std::abs(s1) ) > eps);

        aa += 2 * y;
        dgauss += s1;
      }
    return dgauss;
  }

  //_________________________________________________________________________
  double Integrator::integrateGK(double const& a, double const& b, double const& eps) const
  {
    const double delta = eps15 * std::abs(a - b);
    double dgauss = 0;
    double aa = a;

    while (true)
      {
        double y = b - aa;
        if (std::abs(y) < delta)
          break;

        double s0;
        double s1;
        do
          {
            const double bb = aa + y;

            // Central point and half interval
            const double c1 = ( bb + aa ) / 2;
            const double c2 = ( bb - aa ) / 2;

            const double f = _func(c1);
            s0 = gk_w[0][0] * f;
            s1 = gk_w[1][0] * f;
            for (int i = 1; i < 4; i++)
              {
                // Integral with smaller number of points
                const double u0 = gk_x[0][i] * c2;
                const double fu = _func(c1 + u0);
                const double fd = _func(c1 - u0);
                s0 += gk_w[0][i] * ( fu + fd );

                // Integral with larger number of points (exploit the
                // function calls used for s0).
                const int j = 2 * i - 1;
                const double u1 = gk_x[1][j] * c2;
                s1 += gk_w[1][j] * ( _func(c1 + u1) + _func(c1 - u1) ) + gk_w[1][j+1] * ( fu + fd );
              }
            double w = gk_x[1][7] * c2;
            s1 += gk_w[1][7] * ( _func(c1 + w) + _func(c1 - w) );
            s0 *= c2;
            s1 *= c2;

            // Split the interval into two
            y /= 2;

            // Stop if the interval is smaller than the minimum step
            if (std::abs(y) < delta)
              throw std::runtime_error(error("Integrator::integrateGK", "Too high accuracy required."));

            // Check accuracy and continue if necessary
          }
        while (std::abs(s1 - s0) / ( 1 + std::abs(s1) ) > eps);

        aa += 2 * y;
        dgauss += s1;
      }
    return dgauss;
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, std::vector<double> const& FixPts, double const& eps) const
  {
    // Create vector of fixed points including the integration bounds
    // in ordered sequence and removing the points outside the
    // integration range.
    const double ap = (a < b ? a : b);
    const double bp = (a < b ? b : a);
    std::vector<double> bounds{ap, bp};
    for (auto const& ib : FixPts)
      if(ib > ap && ib < bp)
        bounds.push_back(ib);

    // Sort vector according on how the integration bounds are
    // ordered.
    std::sort(bounds.begin(), bounds.end(), [a, b] (double const& x1, double const& x2) -> bool { return (b > a ? x1 < x2 : x1 > x2); });

    // Now compute sub-integrals and sum them up
    double dgauss = 0;
    for (int i = 1; i < (int) bounds.size(); i++)
      dgauss += integrate(bounds[i-1], bounds[i], eps);

    return dgauss;
  }
}
