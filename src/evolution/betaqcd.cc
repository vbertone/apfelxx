//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/betaqcd.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double beta0(int const& nf)
  {
    const double coeff = ( 11 * CA - 4 * TR * nf ) / 3;
    return coeff;
  }

  //_________________________________________________________________________
  double beta1(int const& nf)
  {
    const double coeff = 34 * CA * CA / 3 - 20 * CA * TR * nf / 3 - 4 * CF * TR * nf;
    return coeff;
  }

  //_________________________________________________________________________
  double beta2(int const& nf)
  {
    const double coeff = 2857 * pow(CA,3) / 54
      + ( 2 * CF * CF - 205 * CF * CA / 9 - 1415 * CA * CA / 27 ) * TR * nf
      + ( 44 * CF / 9 + 158 * CA / 27 ) * TR * TR * nf * nf;
    return coeff;
  }

  //_________________________________________________________________________
  double beta3(int const& nf)
  {
    const double coeff = 149753. / 6. + 3564. * zeta3
      + ( - 1078361. / 162. - 6508. * zeta3 / 27. ) * nf
      + ( 50065. / 162. + 6472. * zeta3 / 81 ) * nf * nf
      + 1093. / 729. * nf * nf * nf;
    return coeff;
  }
}
