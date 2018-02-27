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
   * @name Special functions
   * Collection of special functions needed in the evaluation of some
   * expressions.
   */
  ///@{
  /**
   * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
   * @param x: the real argument
   * @return \f$\mathrm{Li}_2(z)\f$
   * @note Implementation translated by R.Brun from CERNLIB DILOG function C332.
   */
  double dilog(double const& x);

  /**
   * @brief function for the computation of the Nielsen's generalized dilogs.
   * @param n: integer argument
   * @param p: integer argument
   * @param x: real argument
   * @return \f$\mathrm{S}_{n,p}(x)\f$
   * @note Implementation translated from CERNLIB WGPLG.
   */
  double wgplg(int const& n, int const& p, double const& x);
  ///@}
}
