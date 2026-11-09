//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscff.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  ReCFF20ns::ReCFF20ns():
    Expression{}
  {
  }
  double ReCFF20ns::Regular(double const& x)  const
  {
    return - 2 / ( 1 + x );
  }
  double ReCFF20ns::Singular(double const& x) const
  {
    return - 2 / ( 1 - x );
  }
  double ReCFF20ns::LocalPP(double const& c, double const& d) const
  {
    return - 2 * ( log(1 - c) - (d > 1 ? log(d - 1) : 0) );
  }
}
