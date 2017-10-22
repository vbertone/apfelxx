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
   * @brief Return the number of active flavours at the scale Q given
   * the (ordered) vector of thresholds.
   * @param Q scale
   * @param Thresholds
   * @return number of active flavours
   */
  int NF(double const& Q, vector<double> const& Thresholds);
}
