//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <functional>
using std::function;

namespace apfel
{
  /**
   * @brief rk4 using lambdas, this solves:
   * y' = f(t,y)
   * where
   * dy = rk4(f(t,y))
   * so differentiation between lower and upper
   * y += dy(k,y,dlk) for k in steps
   */
  function<double (double const&, double const&, double const&)>
  rk4(function<double(double const&,double const&)> const& f);

}
