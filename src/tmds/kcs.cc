//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/kcs.h"
#include "apfel/betaqcd.h"
#include "apfel/gammak.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double KCS00()
  {
    return 0;
  }

  //_________________________________________________________________________
  double KCS01()
  {
    return - gammaK0();
  }

  //_________________________________________________________________________
  double KCS10(int const& nf)
  {
    return CA * ( 28 * zeta3 - 808. / 27. ) + 224 * TR * nf / 27;
  }

  //_________________________________________________________________________
  double KCS11(int const& nf)
  {
    return - gammaK1(nf);
  }

  //_________________________________________________________________________
  double KCS12(int const& nf)
  {
    return beta0qcd(nf) * gammaK0();
  }

  //_________________________________________________________________________
  double KCS20(int const& nf)
  {
    return 2 * ( CA * CA * ( - 176 * zeta3 * zeta2 / 3 + 6392 * zeta2 / 81 + 12328 * zeta3 / 27
                             + 154 * zeta4 / 3 - 192 * zeta5 - 297029. / 729. ) / 2
                 + CA * TR * nf * ( - 824 * zeta2 / 81 - 904 * zeta3 / 27 + 20 * zeta4 / 3 + 62626. / 729. )
                 + 2 * TR * TR * nf * nf * ( - 32 * zeta3 / 9 - 1856. / 729. )
                 + CF * TR * nf * ( - 304 * zeta3 / 9 - 16 * zeta4 + 1711. / 27. ) );
  }

  //_________________________________________________________________________
  double KCS21(int const& nf)
  {
    return - gammaK2(nf);
  }

  //_________________________________________________________________________
  double KCS22(int const& nf)
  {
    return beta1qcd(nf) * gammaK0() + 2 * beta0qcd(nf) * gammaK1(nf);
  }

  //_________________________________________________________________________
  double KCS23(int const& nf)
  {
    const double b02 = pow(beta0qcd(nf), 2);
    return - 4 * b02 * gammaK0() / 3;
  }
}
