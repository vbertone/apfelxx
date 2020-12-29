//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <functional>
#include <vector>

namespace apfel
{
  /**
   * @brief The Integrator class performs unidimensional numerical
   * integrations using the Guassian quadrature.
   */
  class Integrator
  {
  public:
    /**
     * @name Enumerator for the currently available integration
     * methods.
     */
    enum IntegrationMethod: int {GAUSS_LEGENDRE, GAUSS_KRONROD};

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The Integrator constructor.
     * @param func: The function of one variable to be integrated
     * @param method: The integration method to be used (default: GAUSS_KRONROD)
     */
    Integrator(std::function<double(double const&)> const& func, IntegrationMethod const& method = GAUSS_KRONROD);

    /**
     * @brief Function that integrates the integrand with a given
     * relative accuracy using the method defined in the constructor.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @param eps: the required relative accuracy
     * @return the value of the integral
     */
    double integrate(double const& xmin, double const& xmax, double const& eps) const;

    /**
     * @brief Function that integrates the integrand using the Gauss-Legendre method.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @return a pair containing the value of the integral computed
     * with the 16-point method and the relative difference w.r.t. the
     * 8-point one.
     */
    std::pair<double, double> integrateGL(double const& xmin, double const& xmax) const;

    /**
     * @brief Function that integrates the integrand using the Gauss-Kronrod method.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @return a pair containing the value of the integral computed
     * with the 15-point method and the relative difference w.r.t. the
     * 7-point one.
     */
    std::pair<double, double> integrateGK(double const& xmin, double const& xmax) const;

    /**
     * @brief Function that integrates the integrand with a given
     * relative accuracy using a set of fixed point on the integration
     * range.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @param FixPts: the vector of fixed points of the integration
     * @param eps: the required relative accuracy
     * @return the value of the integral
     */
    double integrate(double const& xmin, double const& xmax, std::vector<double> const& FixPts, double const& eps) const;

    /**
     * @brief Function for the integrand.
     * @param x: the integration variable
     * @return the integrand evaluated at x
     */
    double integrand(double const& x) const { return _func(x); };

    /**
     * @brief Function that returns the integration method.
     */
    IntegrationMethod Method() const { return _method; };

  private:
    std::function<double(double const&)> const _func;   //!< The integrand function
    IntegrationMethod                    const _method; //!< The integration method
  };
}
