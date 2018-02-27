//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>

using std::vector;

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
  int NF(double const& Q, vector<double> const& Thresholds);
  ///@}
}
