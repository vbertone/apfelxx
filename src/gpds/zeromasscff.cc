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
  ReCFF10ns::ReCFF10ns():
    Expression{}
  {
  }
  double ReCFF10ns::Regular(double const& x)  const
  {
    return - 1 / ( 1 + x );
  }
  double ReCFF10ns::Singular(double const& x) const
  {
    return - 1 / ( 1 - x );
  }
  double ReCFF10ns::LocalPP(double const& x)    const
  {
    return - log(1 - x);
  }

  //_________________________________________________________________________________
  ImCFF11ns::ImCFF11ns():
    Expression{}
  {
  }

  //_________________________________________________________________________________
  ReCFF11ns::ReCFF11ns():
    Expression{}
  {
  }

  //_________________________________________________________________________________
  ImCFF11g::ImCFF11g():
    Expression{}
  {
  }

  //_________________________________________________________________________________
  ReCFF11g::ReCFF11g():
    Expression{}
  {
  }
}
