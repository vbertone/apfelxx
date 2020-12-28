//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <functional>
#include <vector>
#include <string>

namespace apfel
{
  /**
   * @brief The Integrator2D class performs two-dimensional numerical
   * integrations using the Guassian quadrature.
   */
  class Integrator2D
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The Integrator constructor.
     * @param func: The function of two variables to be integrated
     */
    Integrator2D(std::function<double(double const&, double const&)> const& func);

    /**
     * @brief Function that integrates the integrand with a given
     * relative accuracy using the method defined in the constructor.
     * @param xmin: the lower bound integration bound for the first variable
     * @param xmax: the upper bound integration bound for the first variable
     * @param ymin: the lower bound integration bound for the second variable
     * @param ymax: the upper bound integration bound for the second variable
     * @param eps: the required relative accuracy
     * @return the value of the integral
     */
    double integrate(double const& xmin, double const& xmax, double const& ymin, double const& ymax, double const& eps) const;

    /**
     * @brief Function that integrates the integrand with a using the
     * 8x8 and 16x16 Gauss-Legendre method.
     * @param xmin: the lower bound integration bound for the first variable
     * @param xmax: the upper bound integration bound for the first variable
     * @param ymin: the lower bound integration bound for the second variable
     * @param ymax: the upper bound integration bound for the second variable
     * @return a pair containing the value of the integral computed
     * with the 16x16 method and the relative difference w.r.t. the
     * 8x8 one.
     */
    std::pair<double, double> integrate(double const& xmin, double const& xmax, double const& ymin, double const& ymax) const;

    /**
     * @brief Function for the integrand.
     * @param x: the first variable
     * @param y: the second variable
     * @return the integrand evaluated at x and y
     */
    double integrand(double const& x, double const& y) const { return _func(x, y); };

  private:
    std::function<double(double const&, double const&)> _func;   //!< The integrand function
  };
}
