//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/distribution.h"

#include <functional>
#include <map>
#include <map>

using namespace std;

namespace apfel
{
  /**
   * @brief Collection of funcitions to rotate distributions from one
   * basis to the other.
   */
  map<int,double> PhysToQCDEv(map<int,double> const& InPhysMap);
  map<int,double> QCDEvToPhys(map<int,double> const& QCDEvMap);

  map<int,Distribution> QCDEvToPhys(map<int,Distribution> const& QCDEvMap);
}
