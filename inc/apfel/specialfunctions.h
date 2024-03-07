//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <map>

namespace apfel
{
  /**
   * @name Fortran harmonic polylogarithms
   * @brief Harmonic polylogarithms up to weight five
   * @param x: real input argument
   * @param nw: maximum number of weights requested
   * @param Hr1: weight 1 harmonic polylogs (1D array)
   * @param Hr2: weight 2 harmonic polylogs (2D array)
   * @param Hr3: weight 3 harmonic polylogs (3D array)
   * @param Hr4: weight 4 harmonic polylogs (4D array)
   * @param Hr5: weight 5 harmonic polylogs (5D array)
   * @param n1: lower bound of the weight index requested
   * @param n2: upper bound of the weight index requested
   * @note This is just a suitably formatted wrapper of the original
   * fortran function (see src/kernel/hplog.f) to facilitate the call
   * of the harmonic logarithms from a C++ code.
   */
  extern "C"
  {
    double apf_hplog_(double *wx, int *wnw, double *Hr1, double *Hr2, double *Hr3, double *Hr4, double *Hr5, int *wn1, int *wn2);
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
   * @brief Function for the computation of the Nielsen's generalized dilogs.
   * @param n: integer argument
   * @param p: integer argument
   * @param x: real argument
   * @return \f$\mathrm{S}_{n,p}(x)\f$
   * @note Implementation translated from CERNLIB WGPLG.
   */
  double wgplg(int const& n, int const& p, double const& x);

  /**
   * @brief Function for the computation of the Harmonic polylogs up
   * to weight 5.
   * @param w: vector of weights
   * @param x: real argument
   * @return \f$\mathrm{H}(\{w\},x)\f$
   * @note C++ adaptation of the FORTRAN implementation discussed in
   * https://arxiv.org/pdf/1809.07084.pdf. The argument x is limited
   * to the interval [0, sqrt(2)-1].
   */
  double hpoly(std::vector<int> const& w, double const& x);

  /**
   * @brief Digamma function.
   * @param x: real argument
   * @return The digamma funciton computed at x.
   * @note C++ (real) adaptation of the FORTRAN (complex)
   * implementation present in the PEGASUS code (hep-ph/0408244).
   */
  double digamma(double const& x);

  /**
   * @brief Function that returns the index to be used with
   * unidimensional arrays returned by hplog_.
   * @param w: the packed vector of weights
   */
  int HPLogMap(std::vector<int> const& w);

  /**
   * @brief Function that returns the unpacked weights of the HPL
   * given the input vector.
   * @param w: the packed vector of weights
   */
  std::vector<int> UnpackWeights(std::vector<int> const& w);
  ///@}
}
