//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/hardfactors.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double H1DY()
  {
    return 2 * CF * ( - 8 + 7 * Pi2 / 6 );
  }

  //_________________________________________________________________________
  double H2DY(int const& nf)
  {
    return 2 * CF * ( CF * ( 511. / 8. - 83 * Pi2 / 6 - 30 * zeta3 + 67 * Pi2 * Pi2 / 60 ) +
                      CA * ( - 51157. / 648. + 1061 * Pi2 / 108 + 313 * zeta3 / 9 - 4 * Pi2 * Pi2 / 45 ) +
                      TR * nf * ( 4085. / 162. - 91 * Pi2 / 27 + 4 * zeta3 / 9 ) );
  }

  //_________________________________________________________________________
  double H1SIDIS()
  {
    return 2 * CF * ( - 8 + zeta2 );
  }

  //_________________________________________________________________________
  double H2SIDIS(int const& nf)
  {
    return 2 * CF * ( CF * ( 511. / 8. + 13 * zeta2 - 30 * zeta3 + 39 * zeta4 / 2 ) +
                      CA * ( - 51157. / 648. - 337 * zeta2 / 18 + 313 * zeta3 / 9 + 22 * zeta4 ) +
                      TR * nf * ( 4085. / 162. + 46 * zeta2 / 9 + 4 * zeta3 / 9 ) );;
  }

  //_________________________________________________________________________
  double H1ggH()
  {
    return 2 * ( CA * ( 5 + 7 * Pi2 / 6 ) - 3 * CF );
  }

  //_________________________________________________________________________
  double H2ggH(int const& nf)
  {
    const double mH = 125.1;
    const double mt = 172.9;
    const double Lt = 2 * log( mH / mt );
    return ( - 135 * CA + 23827 * CA * CA - 216 * CF - 15660 * CA * CF + 5832 * CF * CF
             + 2268 * CA * CA * Lt - 3564 * CA * CF * Lt - 4510 * CA * nf - 4428 * CF * nf
             + 1296 * CF * nf * Lt + 6795 * CA * CA * Pi2 - 2268 * CA * CF * Pi2
             - 450 * Pi2 * CA * nf + 333 * Pi2 * Pi2 * CA * CA - 5148 * CA * CA * zeta3
             - 1656 * zeta3 * CA * nf + 2592 * zeta3 * CF * nf ) / 162;;
  }
}
