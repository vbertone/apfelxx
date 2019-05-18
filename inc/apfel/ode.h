//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Runge-Kutta (RK) ODE solvers.
   * These functions solve the ordinary differential equation (ODE):
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
   */
  ///@{
  /**
   * @brief Template function that implements the fourth order RK
   * algorithm.
   * @param f: the function on the r.h.s. of the ODE
   * @return the function tha returns the step
   */
  template<class U>
  std::function<U(double const&, U const&, double const&)>
  rk4(std::function<U(double const& t, U const& Obj)> const& f)
  {
    // *INDENT-OFF*
    return
      [       f            ](double const& t, U const& y,  double const& dt) -> U{ return
      [t,y,dt,f            ](                 U const& dy1                 ) -> U{ return
      [t,y,dt,f,dy1        ](                 U const& dy2                 ) -> U{ return
      [t,y,dt,f,dy1,dy2    ](                 U const& dy3                 ) -> U{ return
      [       f,dy1,dy2,dy3](                 U const& dy4                 ) -> U{ return
      ( dy1 + 2 * dy2 + 2 * dy3 + dy4 ) / 6  ;} (
      dt * f( t + dt    , y + dy3     )     );} (
      dt * f( t + dt / 2, y + dy2 / 2 )     );} (
      dt * f( t + dt / 2, y + dy1 / 2 )     );} (
      dt * f( t         , y           )     );} ;
    // *INDENT-ON*
  }

  /**
   * @brief Template function that implements the first order RK
   * algorithm.
   * @param f: the function on the r.h.s. of the ODE
   * @return the function tha returns the step
   */
  template<class U>
  std::function<U(double const&, U const&, double const&)>
  rk1(std::function<U(double const& t, U const& Obj)> const& f)
  {
    return [f](double const& t, U const& y,  double const& dt) -> U{ return dt * f(t, y); } ;
  }
  ///@}
}
