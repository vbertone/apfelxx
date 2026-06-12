//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
    return [f] (double const& t, U const& y, double const& dt) -> U
    {
      const U dy1 = dt * f( t, y           );
      const U dy2 = dt * f( t + dt / 2, y + dy1 / 2 );
      const U dy3 = dt * f( t + dt / 2, y + dy2 / 2 );
      const U dy4 = dt * f( t + dt, y + dy3     );
      return ( dy1 + 2 * dy2 + 2 * dy3 + dy4 ) / 6;
    };
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
    return [f] (double const& t, U const& y,  double const& dt) -> U{ return dt * f(t, y); } ;
  }
  ///@}
}
