//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/set.h"
#include "apfel/distribution.h"

#include <functional>

using std::function;

namespace apfel
{
  /**
   * @brief rk4 using lambdas, this solves:
   *
   *    dy / dt = f(t,y)
   *
   * where:
   *
   *    dy = rk4(f(t,y))
   *
   * so differentiation between lower and upper:
   *
   *    y += dy(t,y,dt)
   *
   * U is the type of the 'y' object.
   * T is the f() signature auto-deduced.
   */
  template<class U = double, class T>
  function<U(double const&, U const&, double const&)>
  rk4(T const& f)
  {
    return
      [       f            ](double const& t, U const& y,  double const& dt) -> U{ return
      [t,y,dt,f            ](                 U const& dy1                 ) -> U{ return
      [t,y,dt,f,dy1        ](                 U const& dy2                 ) -> U{ return
      [t,y,dt,f,dy1,dy2    ](                 U const& dy3                 ) -> U{ return
      [t,y,dt,f,dy1,dy2,dy3](                 U const& dy4                 ) -> U{ return
      ( dy1 + 2 * dy2 + 2 * dy3 + dy4 ) / 6  ;} (
      dt * f( t + dt    , y + dy3     )     );} (
      dt * f( t + dt / 2, y + dy2 / 2 )     );} (
      dt * f( t + dt / 2, y + dy1 / 2 )     );} (
      dt * f( t         , y           )     );} ;
  }
}
