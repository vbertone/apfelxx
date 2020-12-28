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
  Integrator2D::Integrator2D(std::function<double(double const&, double const&)> const& func, Integrator::IntegrationMethod const& method):
    _func(func),
    _method(method)
  {
    switch (_method)
      {
      case Integrator::GAUSS_LEGENDRE:
        _integrate = std::bind(&Integrator2D::integrateGL, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        break;
      case Integrator::GAUSS_KRONROD:
        _integrate = std::bind(&Integrator2D::integrateGK, this, std::placeholders::_1, std::placeholders::_2, std::placeholders::_3, std::placeholders::_4);
        break;
      default:
        throw std::runtime_error(error("Integrator2D::Integrator2D", "Unknown integration method"));
      }
  }

  //_________________________________________________________________________
  double Integrator2D::integrate(double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& eps) const
  {
    // Stop if the interval is too small
    if (std::abs( ( xmax - xmin ) / ( xmax + xmin ) ) < eps15 || std::abs( ( ymax - ymin ) / ( ymax + ymin ) ) < eps15)
      throw std::runtime_error(error("Integrator2D::integrate", "Too high accuracy required."));

    const std::pair<double, double> integ = _integrate(xmin, xmax, ymin, ymax);
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

  //_________________________________________________________________________
  std::pair<double, double> Integrator2D::integrateGL(double const& a, double const& b, double const& c, double const& d) const
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
  std::pair<double, double> Integrator2D::integrateGK(double const& a, double const& b, double const& c, double const& d) const
  {
    // Central point and half interval
    const double c1x = ( b + a ) / 2;
    const double c2x = ( b - a ) / 2;
    const double c1y = ( d + c ) / 2;
    const double c2y = ( d - c ) / 2;

    double const f00 = _func(c1x, c1y);
    double s0 = gk_w[0][0] * gk_w[0][0] * f00;
    double s1 = gk_w[1][0] * gk_w[1][0] * f00;
    for (int j = 1; j < 4; j++)
      {
        const double fu0 = _func(c1x + gk_x[0][j] * c2x, c1y);
        const double fd0 = _func(c1x - gk_x[0][j] * c2x, c1y);
        const double f0u = _func(c1x, c1y + gk_x[0][j] * c2y);
        const double f0d = _func(c1x, c1y - gk_x[0][j] * c2y);
        s0 += gk_w[0][0] * gk_w[0][j] * ( fu0 + fd0 + f0u + f0d );

        const int k = 2 * j - 1;
        s1 += gk_w[1][0] * gk_w[1][k] * ( _func(c1x + gk_x[1][k] * c2x, c1y) + _func(c1x - gk_x[1][k] * c2x, c1y) +
                                          _func(c1x, c1y + gk_x[1][k] * c2y) + _func(c1x, c1y - gk_x[1][k] * c2y) )
              + gk_w[1][0] * gk_w[1][k+1] * ( fu0 + fd0 + f0u + f0d );
      }
    s1 += gk_w[1][0] * gk_w[1][7] * ( _func(c1x + gk_x[1][7] * c2x, c1y) + _func(c1x - gk_x[1][7] * c2x, c1y) +
                                      _func(c1x, c1y + gk_x[1][7] * c2y) + _func(c1x, c1y - gk_x[1][7] * c2y) );

    for (int i = 1; i < 4; i++)
      {
        const int l = 2 * i - 1;
        const double ux1 = gk_x[1][l] * c2x;
        const double ux2 = gk_x[1][l+1] * c2x;
        for (int j = 1; j < 4; j++)
          {
            const int k = 2 * j - 1;
            const double uy1 = gk_x[1][k] * c2y;
            const double uy2 = gk_x[1][k+1] * c2y;
            const double fuu = _func(c1x + ux2, c1y + uy2);
            const double fud = _func(c1x + ux2, c1y - uy2);
            const double fdu = _func(c1x - ux2, c1y + uy2);
            const double fdd = _func(c1x - ux2, c1y - uy2);
            s0 += gk_w[0][i] * gk_w[0][j] * ( fuu + fud + fdu + fdd );

            s1 += gk_w[1][l] * gk_w[1][k] * ( _func(c1x + ux1, c1y + uy1) + _func(c1x + ux1, c1y - uy1) +
                                              _func(c1x - ux1, c1y + uy1) + _func(c1x - ux1, c1y - uy1) )
                  + gk_w[1][l] * gk_w[1][k+1] * ( _func(c1x + ux1, c1y + uy2) + _func(c1x + ux1, c1y - uy2) +
                                                  _func(c1x - ux1, c1y + uy2) + _func(c1x - ux1, c1y - uy2) )
                  + gk_w[1][l+1] * gk_w[1][k] * ( _func(c1x + ux2, c1y + uy1) + _func(c1x + ux2, c1y - uy1) +
                                                  _func(c1x - ux2, c1y + uy1) + _func(c1x - ux2, c1y - uy1) )
                  + gk_w[1][l+1] * gk_w[1][k+1] * ( fuu + fud + fdu + fdd );
          }

        const double uy1 = gk_x[1][7] * c2y;
        const double uy2 = gk_x[1][7] * c2y;
        s1 += gk_w[1][l] * gk_w[1][7] * ( _func(c1x + ux1, c1y + uy1) + _func(c1x + ux1, c1y - uy1) +
                                          _func(c1x - ux1, c1y + uy1) + _func(c1x - ux1, c1y - uy1) )
              + gk_w[1][l+1] * gk_w[1][7] * ( _func(c1x + ux2, c1y + uy2) + _func(c1x + ux2, c1y - uy2) +
                                              _func(c1x - ux2, c1y + uy2) + _func(c1x - ux2, c1y - uy2) );
      }
    const double ux = gk_x[1][7] * c2x;
    for (int j = 1; j < 4; j++)
      {
        const int k = 2 * j - 1;
        const double uy1 = gk_x[1][k] * c2y;
        const double uy2 = gk_x[1][k+1] * c2y;
        s1 += gk_w[1][7] * gk_w[1][k] * ( _func(c1x + ux, c1y + uy1) + _func(c1x + ux, c1y - uy1) +
                                          _func(c1x - ux, c1y + uy1) + _func(c1x - ux, c1y - uy1) )
              + gk_w[1][7] * gk_w[1][k+1] * ( _func(c1x + ux, c1y + uy2) + _func(c1x + ux, c1y - uy2) +
                                              _func(c1x - ux, c1y + uy2) + _func(c1x - ux, c1y - uy2) );
      }
    const double uy = gk_x[1][7] * c2y;
    s1 += gk_w[1][7] * gk_w[1][7] * ( _func(c1x + ux, c1y + uy) + _func(c1x + ux, c1y - uy) +
                                      _func(c1x - ux, c1y + uy) + _func(c1x - ux, c1y - uy) );
    s0 *= c2x * c2y;
    s1 *= c2x * c2y;
    return {s1, std::abs(s1 - s0) / ( 1 + std::abs(s1) )};
  }
}
