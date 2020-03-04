//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <functional>
#include <vector>

namespace apfel
{
  /**
   * @brief The Integrator class perform unidimensional numerical
   * integration using the Guassian quadrature.
   */
  class Integrator
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The Integrator constructor.
     * @param func: The function of one variable to be integrated
     */
    Integrator(std::function<double(double const&)> const& func);

    /**
     * @brief The Integrator constructor.
     * @param func: the function of two variables to be integrated over the first
     * @param arg2: the value of the second variable while integrating over the first
     */
    Integrator(std::function<double(double const&, double const&)> const& func2, double const& arg2);

    /**
     * @brief The Integrator constructor.
     * @param func: the function of three variables to be integrated over the first
     * @param arg2: the value of the second variable while integrating over the first
     * @param arg3: the value of the third variable while integrating over the first
     */
    Integrator(std::function<double(double const&, double const&, double const&)> const& func3, double const& arg2, double const& arg3);
    ///@}

    /**
     * @brief Function that integrates the integrand with a given
     * relative accuracy.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @param eps: the required relative accuracy
     * @return the value of the integral
     */
    double integrate(double const& xmin, double const& xmax, double const& eps) const;

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
     * @brief Function that integrates the integrand using a given
     * number of point for the gauss quadrature.
     * @param xmin: the lower bound integration bound
     * @param xmax: the upper bound integration bound
     * @param m: number of point of the Guass quadrature
     * @return the value of the integral
     */
    double integrate(double const& xmin, double const& xmax, int const& m) const;

    /**
     * @brief Function for the integrand.
     * @param x: the integration variable
     * @return the integrand evaluated at x
     */
    double integrand(double const& x) const { return _func(x); };

  private:
    std::function<double(double const&)> _func; //!< The integrand function
  };
}
