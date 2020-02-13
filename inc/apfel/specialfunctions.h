//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Fortran harmonic polylogarithms
   * @brief Harmonic polylogarithms up to weight four
   * @param x: real input argument
   * @param nw: maximum number of weights requested
   * @param Hr1: weight 1 harmonic polylogs (1D array)
   * @param Hr2: weight 2 harmonic polylogs (2D array)
   * @param Hr3: weight 3 harmonic polylogs (3D array)
   * @param Hr4: weight 4 harmonic polylogs (4D array)
   * @param n1: lower bound of the weight index requested
   * @param n2: upper bound of the weight index requested
   * @note This is just a suitably formatted wrapper of the original
   * fortran function (see src/kernel/hplog.f) to facilitate the call
   * of the harmonic logarithms from a C++ code.
   */
  extern"C"
  {
    double hplog_(double *wx, int *wnw, double *Hr1, double *Hr2, double *Hr3, double *Hr4, int *wn1, int *wn2);
  }

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
