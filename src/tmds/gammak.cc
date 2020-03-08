//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gammak.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double gammaK0()
  {
    return 8;
  }

  //_________________________________________________________________________
  double gammaK1(int const& nf)
  {
    return 8 * ( ( 67. / 9. - Pi2 / 3 ) * CA - 20 * TR * nf / 9 );
  }

  //_________________________________________________________________________
  double gammaK2(int const& nf)
  {
    return 8 * ( ( 245. / 6. - 134 * Pi2 / 27 + 11 * Pi2 * Pi2 / 45 + 22 * zeta3 / 3 ) * CA * CA
                 + ( - 418. / 27. + 40 * Pi2 / 27 - 56 * zeta3 / 3 ) * CA * TR * nf
                 + ( - 55. / 3. + 16 * zeta3 ) * CF * TR * nf
                 - 16 * TR * TR * nf * nf / 27 );
  }

  //_________________________________________________________________________
  double gammaK3(int const& nf)
  {
    const int nf2 = nf * nf;
    const int nf3 = nf2 * nf;
    return ( 20702 - 5171.9 * nf + 195.5772 * nf2 + 3.272344 * nf3 ) / CF;
  }

  //_________________________________________________________________________
  double gammaK3gmq()
  {
    return 40880 / CA - 20702 / CF;
  }
}
