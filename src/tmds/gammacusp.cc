//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/gammacusp.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double GammaCusp0()
  {
    return 1;
  }

  //_________________________________________________________________________
  double GammaCusp1(int const& nf)
  {
    const double coeff = ( 67. / 9. - Pi2 / 3 ) * CA - 20 * TR * nf / 9;
    return coeff;
  }

  //_________________________________________________________________________
  double GammaCusp2(int const& nf)
  {
    const double coeff = ( 245. / 6. - 134 * Pi2 / 27 + 11 * Pi2 * Pi2 / 45 + 22 * zeta3 / 3 ) * CA * CA
      + ( - 418. / 27. + 40 * Pi2 / 27 - 56 * zeta3 / 3 ) * CA * TR * nf
      + ( - 55. / 3. + 16 * zeta3 ) * CF * TR * nf
      - 16 * TR * TR * nf * nf / 27;
    return coeff;
  }
}
