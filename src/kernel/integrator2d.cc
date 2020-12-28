//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/integrator2d.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

#include <iomanip>
#include <map>

namespace apfel
{
  //_________________________________________________________________________
  Integrator2D::Integrator2D(std::function<double(double const&, double const&)> const& func):
    _func(func)
  {
  }

  //_________________________________________________________________________
  std::pair<double, double> Integrator2D::integrate(double const& a, double const& b, double const& c, double const& d) const
  {
    // Central point and half interval
    const double c1x = ( b + a ) / 2;
    const double c2x = ( b - a ) / 2;
    const double c1y = ( d + c ) / 2;
    const double c2y = ( d - c ) / 2;

    // Integral with smaller number of points
    double s0 = 0;
    for (int i = 0; i < (int) gl_x[0].size(); i++)
      {
        const double ux = gl_x[0][i] * c2x;
        double sx = 0;
        for (int j = 0; j < (int) gl_x[0].size(); j++)
          {
            const double uy = gl_x[0][j] * c2y;
            sx += gl_w[0][j] * ( _func(c1x + ux, c1y + uy) +
                                 _func(c1x + ux, c1y - uy) +
                                 _func(c1x - ux, c1y + uy) +
                                 _func(c1x - ux, c1y - uy) );
          }
        s0 += gl_w[0][i] * sx;
      }
    s0 *= c2x * c2y;

    // Integral with larger number of points
    double s1 = 0;
    for (int i = 0; i < (int) gl_x[1].size(); i++)
      {
        const double ux = gl_x[1][i] * c2x;
        double sx = 0;
        for (int j = 0; j < (int) gl_x[1].size(); j++)
          {
            const double uy = gl_x[1][j] * c2y;
            sx += gl_w[1][j] * ( _func(c1x + ux, c1y + uy) +
                                 _func(c1x + ux, c1y - uy) +
                                 _func(c1x - ux, c1y + uy) +
                                 _func(c1x - ux, c1y - uy) );
          }
        s1 += gl_w[1][i] * sx;
      }
    s1 *= c2x * c2y;

    return {s1, std::abs(s1 - s0) / ( 1 + std::abs(s1) )};
  }

  //_________________________________________________________________________
  double Integrator2D::integrate(double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& eps) const
  {
    // Stop if the interval is too small
    if (std::abs( ( xmax - xmin ) / ( xmax + xmin ) ) < eps15 || std::abs( ( ymax - ymin ) / ( ymax + ymin ) ) < eps15)
      throw std::runtime_error(error("Integrator2D::integrate", "Too high accuracy required."));

    const std::pair<double, double> integ = integrate(xmin, xmax, ymin, ymax);
    if (integ.second < eps)
      return integ.first;
    else
      {
        const double xm = ( xmin + xmax ) / 2;
        const double ym = ( ymin + ymax ) / 2;
        return
          integrate(xmin, xm, ymin, ym, eps) +
          integrate(xm, xmax, ymin, ym, eps) +
          integrate(xmin, xm, ym, ymax, eps) +
          integrate(xm, xmax, ym, ymax, eps);
      }
  }
}
