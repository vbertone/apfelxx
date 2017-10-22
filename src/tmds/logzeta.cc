//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/logzeta.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double Lzetaq10()
  {
    return - 3. / 2.;
  }

  //_________________________________________________________________________
  double Lzetaq11()
  {
    return 1. / 2.;
  }

  //_________________________________________________________________________
  double Lzetaq20(int const& nf)
  {
    const double coeff = CF * ( - 3. / 4. + Pi2 - 12 * zeta3 )
      + CA * ( 649. / 108. - 17 * Pi2 / 12 + 19 * zeta3 / 2 )
      + TR * nf * ( - 53. / 27. + Pi2 / 3 );
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetaq22(int const& nf)
  {
    const double coeff = ( 11 * CA - 4 * TR * nf ) / 36;
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetag10(int const& nf)
  {
    const double coeff = - 11. / 6. + 2 * TR * nf / 3 / CA;
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetag11()
  {
    return 1. / 2.;
  }

  //_________________________________________________________________________
  double Lzetag20(int const& nf)
  {
    const double coeff = CA * ( 247. / 54. - 11 * Pi2 / 36 - 5 * zeta3 / 2 )
      + TR * nf * ( - 16. / 3. + Pi2 / 9 )
      + TR * nf * ( 2 * CF + 40 * TR * nf / 27 ) / CA;
    return coeff;
  }

  //_________________________________________________________________________
  double Lzetag22(int const& nf)
  {
    const double coeff = ( 11 * CA - 4 * TR * nf ) / 36;
    return coeff;
  }
}
