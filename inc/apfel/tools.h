//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>

namespace apfel
{
  /**
   * @name Tools
   * Collection of useful tools.
   */
  ///@{
  /**
   * @brief Return the number of active flavours at the scale Q given
   * the (ordered) vector of thresholds.
   * @param Q: the scale
   * @param Thresholds: the vector of thresholds
   * @return number of active flavours at Q
   */
  int NF(double const& Q, std::vector<double> const& Thresholds);

  /**
   * @brief Utility function used by the heavy-quark initiated massive
   * coefficient functions.


   * @param a: first parameter
   * @param b: second parameter
   * @param c: thied parameter
   * @return Triangular function
   */
  double DeltaFun(double const& a, double const& b, double const& c);
  ///@}
}
