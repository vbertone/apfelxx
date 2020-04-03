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
    const int nf2       = nf * nf;
    const int nf3       = nf2 * nf;
    const double CF2    = CF * CF;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double zeta32 = zeta3 * zeta3;
    const double dFF4nc = 5. / 36.;
    const double dFA4nc = 5. / 2.;
    return
      -(nf3*(0.3950617283950617 - (64*zeta3)/27.))
      + CA*nf2*(11.395061728395062 - (608*zeta2)/81. + (2240*zeta3)/27. - (112*zeta4)/3.)
      + CF*nf2*(29.530864197530864 - (640*zeta3)/9. + 32*zeta4)
      + dFF4nc/CF*nf*(256*zeta2 - (256*zeta3)/3. - (1280*zeta5)/3.)
      + CF2*nf*(63.55555555555556 + (592*zeta3)/3. - 320*zeta5)
      + CA*CF*nf*(-420.5679012345679 + (440*zeta2)/3. + (3712*zeta3)/9. - 128*zeta2*zeta3
                  - 176*zeta4 + 160*zeta5)
      + CA2*nf*(-297.98765432098764 + (20320*zeta2)/81. - (23104*zeta3)/27.
                + (448*zeta2*zeta3)/3. - (176*zeta4)/3. + (2096*zeta5)/9.)
      + dFA4nc/CF*(-128*zeta2 + (128*zeta3)/3. - 384*zeta32 + (3520*zeta5)/3. - 992*zeta6)
      + CA3*(1040.469135802469 - (88400*zeta2)/81. + (20944*zeta3)/27. - (352*zeta2*zeta3)/3.
             - 16*zeta32 + 1804*zeta4 - (3608*zeta5)/9. - (2504*zeta6)/3.);
  }

  //_________________________________________________________________________
  double gammaK3gmq()
  {
    return 40880 / CA - 20702 / CF;
  }
}
