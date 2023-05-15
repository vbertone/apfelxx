//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/betaqcd.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double beta0qcd(int const& nf)
  {
    return ( 11 * CA - 4 * TR * nf ) / 3;
  }

  //_________________________________________________________________________
  double beta1qcd(int const& nf)
  {
    return 34 * CA * CA / 3 - 20 * CA * TR * nf / 3 - 4 * CF * TR * nf;
  }

  //_________________________________________________________________________
  double beta2qcd(int const& nf)
  {
    return 2857 * pow(CA,3) / 54
           + ( 2 * CF * CF - 205 * CF * CA / 9 - 1415 * CA * CA / 27 ) * TR * nf
           + ( 44 * CF / 9 + 158 * CA / 27 ) * TR * TR * nf * nf;
  }

  //_________________________________________________________________________
  double beta3qcd(int const& nf)
  {
    return 149753. / 6. + 3564. * zeta3
           + ( - 1078361. / 162. - 6508. * zeta3 / 27. ) * nf
           + ( 50065. / 162. + 6472. * zeta3 / 81 ) * nf * nf
           + 1093. / 729. * nf * nf * nf;
  }

  //_________________________________________________________________________
  double beta4qcd(int const& nf)
  {
    return 8157455. / 16. + 621885. / 2. * zeta3 - 88209. / 2. * zeta4 - 288090. * zeta5
           + ( - 336460813. / 1944. - 4811164. / 81. * zeta3 + 33935. / 6. * zeta4 + 1358995. / 27. * zeta5 ) * nf
           + ( 25960913. / 1944. + 698531. / 81. * zeta3 - 10526. / 9. * zeta4 - 381760. / 81. * zeta5 ) * pow(nf, 2)
           + ( - 630559. / 5832. - 48722. / 243. * zeta3 + 1618. / 27. * zeta4 + 460. / 9. * zeta5 ) * pow(nf, 3)
           + ( 1205. / 2916. - 152 / 81. * zeta3 ) * pow(nf, 4);
  }
}
