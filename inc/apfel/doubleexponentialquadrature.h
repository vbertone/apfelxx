//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <functional>

namespace apfel
{
  /**
   * @brief DE-Quadrature
   * Numerical automatic integrator for improper integral using double
   * dxponential (DE) quadrature. The code is a manipulation of the
   * code linked here:
   *
   * http://www.kurims.kyoto-u.ac.jp/~ooura/intde.html
   */
  class DoubleExponentialQuadrature
  {
  public:
    /**
     * @brief The Integrator constructor.
     * @param eps: relative integration accuracy
     */
    DoubleExponentialQuadrature(double const& eps = 1e-5);

    /**
     * @brief Function that transform the input function.
     * @param f: function to be transformed
     * @param qT: value of qT in which to compute the transform
     * @return the value of the transform
     */
    template<class T>
    T transform(std::function<T(double const&)> const& f, double const& qT) const;

  private:
    double _aw[8000];
  };
}
