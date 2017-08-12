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
   * @brief Real dilogarithm \f$\mathrm{Li}_2(z)\f$
   * @param x real argument
   * @note Implementation translated by R.Brun from CERNLIB DILOG function C332
   * @return \f$\mathrm{Li}_2(z)\f$
   */
  double dilog(double const& x);

  /**
   * @brief function for the computation of the Nielsen's generalized dilogs.
   * @param n integer argument
   * @param p integer argument
   * @param x real argument
   * @note Implementation translated from CERNLIB WGPLG
   * @return \f$\mathrm{S}_{n,p}(x)\f$
   */
  double wgplg(int const& n, int const& p, double const& x);
}
