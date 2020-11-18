//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscoefficientfunctionspol_sl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  G41ns::G41ns():
    Expression()
  {
  }
  double G41ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 3 + 2 * x );
  }
  double G41ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double G41ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x), 2) - 3 * log( 1 - x ) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  GL1ns::GL1ns():
    Expression()
  {
  }
  double GL1ns::Regular(double const& x) const
  {
    return 4 * CF * x;
  }

  //_________________________________________________________________________________
  G11ns::G11ns():
    Expression()
  {
  }
  double G11ns::Regular(double const& x) const
  {
    return 2 * CF * ( - ( 1 + x ) * log( 1 - x ) - ( 1 + pow(x,2) ) * log(x) / ( 1 - x ) + 2 + x );
  }
  double G11ns::Singular(double const& x) const
  {
    return 2 * CF * ( 2 * log( 1 - x ) - 3 / 2. ) / ( 1 - x );
  }
  double G11ns::Local(double const& x) const
  {
    return 2 * CF * ( pow(log(1-x),2) - 3 * log( 1 - x ) / 2 - ( 2 * zeta2 + 9 / 2. ) );
  }

  //_________________________________________________________________________________
  G11g::G11g():
    Expression()
  {
  }
  double G11g::Regular(double const& x) const
  {
    return 4 * TR * ( ( 2 * x - 1 ) * log( ( 1 - x ) / x ) - 4 * x + 3 );
  }
}
