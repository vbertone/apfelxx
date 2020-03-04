//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gammacs.h"
#include "apfel/betaqcd.h"
#include "apfel/gammacusp.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double CSd10()
  {
    return 0;
  }

  //_________________________________________________________________________
  double CSd11()
  {
    return - 4;
  }

  //_________________________________________________________________________
  double CSd20(int const& nf)
  {
    return CA * ( 28 * zeta3 - 808. / 27. ) + 224 * TR * nf / 27;
  }

  //_________________________________________________________________________
  double CSd21(int const& nf)
  {
    return - GammaCusp1(nf);
  }

  //_________________________________________________________________________
  double CSd22(int const& nf)
  {
    return - 2 * beta0qcd(nf);
  }

  //_________________________________________________________________________
  double CSd30(int const& nf)
  {
    return 2 * ( CA * CA * ( - 176 * zeta3 * zeta2 / 3 + 6392 * zeta2 / 81 + 12328 * zeta3 / 27
                             + 154 * zeta4 / 3 - 192 * zeta5 - 297029. / 729. ) / 2
                 + CA * TR * nf * ( - 824 * zeta2 / 81 - 904 * zeta3 / 27 + 20 * zeta4 / 3 + 62626. / 729. )
                 + 2 * TR * TR * nf * nf * ( - 32 * zeta3 / 9 - 1856. / 729. )
                 + CF * TR * nf * ( - 304 * zeta3 / 9 - 16 * zeta4 + 1711. / 27. ) );
  }

  //_________________________________________________________________________
  double CSd31(int const& nf)
  {
    return 2 * beta0qcd(nf) * CSd20(nf) - GammaCusp2(nf);
  }

  //_________________________________________________________________________
  double CSd32(int const& nf)
  {
    return - GammaCusp1(nf) * beta0qcd(nf) - 2 * beta1qcd(nf);
  }

  //_________________________________________________________________________
  double CSd33(int const& nf)
  {
    const double b0 = beta0qcd(nf);
    return - 4 * b0 * b0 / 3;
  }
}
