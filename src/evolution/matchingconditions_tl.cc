//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingconditions_tl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  ATS1Hg_0::ATS1Hg_0():
    Expression()
  {
  }
  double ATS1Hg_0::Regular(double const& x) const
  {
    return 2 * CF * ( 1 + ( 1 - x ) * ( 1 - x ) ) * ( - 1 - 2 * log(x) ) / x;
  }

  //_________________________________________________________________________________
  ATS1Hg_L::ATS1Hg_L():
    Expression()
  {
  }
  double ATS1Hg_L::Regular(double const& x) const
  {
    return 2 * CF * ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
  }

  //_________________________________________________________________________________
  ATS1ggH_L::ATS1ggH_L():
    Expression()
  {
  }
  double ATS1ggH_L::Local(double const&) const
  {
    return - 4 * TR / 3.;
  }

  //_________________________________________________________________________________
  ATS1HH_L::ATS1HH_L():
    Expression()
  {
  }
  double ATS1HH_L::Singular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(x, 2) ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  ATS1HH_0::ATS1HH_0():
    Expression()
  {
  }
  double ATS1HH_0::Singular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(x, 2) ) * ( - 1 - log(1 - x) ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  ATS1gH_L::ATS1gH_L():
    Expression()
  {
  }
  double ATS1gH_L::Regular(double const& x) const
  {
    return 4 * TR * ( x * x + ( 1 - x ) * ( 1 - x ) );
  }
}
