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
  double Integrator::integrate(double const& xmin, double const& xmax, double const& eps) const
  {
    // Stop if the interval is too small
    if (std::abs( ( xmax - xmin ) / ( xmax + xmin ) ) < eps15)
      throw std::runtime_error(error("Integrator::integrate", "Too high accuracy required."));

    // Select integration method
    std::pair<double, double> integ;
    switch (_method)
      {
      case GAUSS_LEGENDRE:
        integ = integrateGL(xmin, xmax);
        break;
      case GAUSS_KRONROD:
        integ = integrateGK(xmin, xmax);
        break;
      default:
        throw std::runtime_error(error("Integrator::integrate", "Unknown integration method"));
      }

    // Recursive call
    if (integ.second < eps)
      return integ.first;
    else
      return integrate(xmin, ( xmin + xmax ) / 2, eps) + integrate(( xmin + xmax ) / 2, xmax, eps);
  }

  //_________________________________________________________________________
  std::pair<double, double> Integrator::integrateGL(double const& a, double const& b) const
  {
    // Central point and half interval
    const double c1 = ( b + a ) / 2;
    const double c2 = ( b - a ) / 2;

    // Integral with smaller number of points
    double s0 = 0;
    for (int i = 0; i < (int) gl_x[0].size(); i++)
      {
        const double u = gl_x[0][i] * c2;
        s0 += gl_w[0][i] * ( _func(c1 + u) + _func(c1 - u) );
      }
    s0 *= c2;

    // Integral with larger number of points
    double s1 = 0;
    for (int i = 0; i < (int) gl_x[1].size(); i++)
      {
        const double u = gl_x[1][i] * c2;
        s1 += gl_w[1][i] * ( _func(c1 + u) + _func(c1 - u) );
      }
    s1 *= c2;

    return {s1, std::abs(s1 - s0) / ( 1 + std::abs(s1) )};
  }

  //_________________________________________________________________________
  std::pair<double, double> Integrator::integrateGK(double const& a, double const& b) const
  {
    // Central point and half interval
    const double c1 = ( b + a ) / 2;
    const double c2 = ( b - a ) / 2;

    const double f = _func(c1);
    double s0 = gk_w[0][0] * f;
    double s1 = gk_w[1][0] * f;
    for (int i = 1; i < (int) gk_x[0].size(); i++)
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

    return {s1, std::abs(s1 - s0) / ( 1 + std::abs(s1) )};
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

  //_________________________________________________________________________
  double Integrator::integrate(std::vector<double> const& FixPts, double const& eps) const
  {
    // Compute sub-integrals and sum them up
    double dgauss = 0;
    for (int i = 1; i < (int) FixPts.size(); i++)
      dgauss += integrate(FixPts[i-1], FixPts[i], eps);
    return dgauss;
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& xmin, double const& xmax, int const& n) const
  {
    // Stop if the index n does not take an allowed value
    if (n != 0 && n != 1)
      throw std::runtime_error(error("Integrator::integrate", "Accuracy index not allowed"));

    // Select integration method
    switch (_method)
      {
      case GAUSS_LEGENDRE:
        return integrateGL(xmin, xmax, n);
      case GAUSS_KRONROD:
        return integrateGK(xmin, xmax, n);
      default:
        throw std::runtime_error(error("Integrator::integrate", "Unknown integration method"));
      }
  }

  //_________________________________________________________________________
  double Integrator::integrateGL(double const& a, double const& b, int const& n) const
  {
    const double c1 = ( b + a ) / 2;
    const double c2 = ( b - a ) / 2;
    double s = 0;
    for (int i = 0; i < (int) gl_x[n].size(); i++)
      {
        const double u = gl_x[n][i] * c2;
        s += gl_w[n][i] * ( _func(c1 + u) + _func(c1 - u) );
      }
    s *= c2;
    return s;
  }

  //_________________________________________________________________________
  double Integrator::integrateGK(double const& a, double const& b, int const& n) const
  {
    const double c1 = ( b + a ) / 2;
    const double c2 = ( b - a ) / 2;
    double s = gk_w[n][0] * _func(c1);
    for (int i = 1; i < (int) gk_x[n].size(); i++)
      {
        const double u = gk_x[n][i] * c2;
        s += gk_w[n][i] * ( _func(c1 + u) + _func(c1 - u) );
      }
    s *= c2;
    return s;
  }

  //_________________________________________________________________________
  double Integrator::integrate(double const& a, double const& b, std::vector<double> const& FixPts, int const& n) const
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
      dgauss += integrate(bounds[i-1], bounds[i], n);

    return dgauss;
  }

  //_________________________________________________________________________
  double Integrator::integrate(std::vector<double> const& FixPts, int const& n) const
  {
    // Compute sub-integrals and sum them up
    double dgauss = 0;
    for (int i = 1; i < (int) FixPts.size(); i++)
      dgauss += integrate(FixPts[i-1], FixPts[i], n);
    return dgauss;
  }
}
