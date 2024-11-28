//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
    return 2 * CF * ( 1 + pow(x, 2) ) * ( - 1 - 2 * log(1 - x) ) / ( 1 - x );
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

  //_________________________________________________________________________________
  ATNS2qqH_0::ATNS2qqH_0():
    Expression()
  {
  }
  double ATNS2qqH_0::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    const double a0   = 4 * ( - 67 * x / 27 + ( - x / 6 + 1. / 3 / ( 1 - x ) - 1. / 6 ) * lnx2
                              + ( - 11 * x / 9 + 10. / 9 / ( 1 - x ) + 1. / 9 ) * lnx + 11. / 27 );
    return CF * TR * a0;
  }
  double ATNS2qqH_0::Singular(double const& x) const
  {
    const double a0 = 224 / 27.;
    return CF * TR * a0 / ( 1 - x );
  }
  double ATNS2qqH_0::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double ln2   = log(2);
    const double ln22  = ln2 * ln2;
    const double ln23  = ln2 * ln22;
    const double a0    = 224 * ln1mx / 27 + 4 * ( 9307. / 648 - 10 * zeta3 / 3 - 19 * Pi2 / 54 - 16 * ln23 / 9
                                                  + 58 * ln22 / 9 + ( 4 * Pi2 / 9 - 359. / 27 ) * ln2 );
    return CF * TR * a0;
  }

  //_________________________________________________________________________________
  ATNS2qqH_L::ATNS2qqH_L():
    Expression()
  {
  }
  double ATNS2qqH_L::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double omeL1 = 8 * ( 1 + x2 ) * lnx / 3 / ( 1 - x ) + 8. / 9 - 88 * x / 9;
    return - CF * TR * omeL1;
  }
  double ATNS2qqH_L::Singular(double const& x) const
  {
    const double omeL1 = 80. / 9 / ( 1 - x );
    return - CF * TR * omeL1;
  }
  double ATNS2qqH_L::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL1 = 80 * ln1mx / 9 + 16 * zeta2 / 3 + 2. / 3;
    return - CF * TR * omeL1;
  }

  //_________________________________________________________________________________
  ATNS2qqH_L2::ATNS2qqH_L2():
    Expression()
  {
  }
  double ATNS2qqH_L2::Regular(double const& x) const
  {
    const double omeL2 = - 4. / 3 - 4 * x / 3;
    return CF * TR * omeL2;
  }
  double ATNS2qqH_L2::Singular(double const& x) const
  {
    const double omeL2 = 8. / 3. / ( 1 - x );
    return CF * TR * omeL2;
  }
  double ATNS2qqH_L2::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL2 = 8 * ln1mx / 3 + 2;
    return CF * TR * omeL2;
  }
}
