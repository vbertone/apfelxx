//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/matchingconditions_tl.h"
#include "apfel/constants.h"

#include <cmath>

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
}
