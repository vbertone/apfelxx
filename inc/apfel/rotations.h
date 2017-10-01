//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

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
  double PhysToQCDEv(int const& i, double const& x, double const& Q, function<double(int const&, double const&, double const&)> const& InDistFunc);
  map<int,double> PhysToQCDEv(double const& x, double const& Q, function<map<int,double>(double const&, double const&)> const& InDistFunc);
  map<int,double> PhysToQCDEv(map<int,double> const& PhysMap);
  map<int,double> QCDEvToPhys(map<int,double> const& QCDEvMap);
}
