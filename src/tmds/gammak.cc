//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
    const double zeta22 = zeta2 * zeta2;
    const double zeta23 = zeta2 * zeta22;
    const double zeta32 = zeta3 * zeta3;
    const double dFF4nc = 5. / 36.;
    const double dFA4nc = 5. / 2.;
    return
      2 * ( nf3 * ( 64. / 27. * zeta3 - 32. / 81. )
            + nf2 * CA * ( - 224. / 15. * zeta22 + 2240. / 27. * zeta3 - 608. / 81. * zeta2 + 923. / 81. )
            + nf2 * CF * ( 64. / 5. * zeta22 - 640. / 9. * zeta3 + 2392. / 81. )
            + nf * CA2 * ( 2096. / 9. * zeta5 + 448. / 3. * zeta3 * zeta2 - 352. / 15. * zeta22
                           - 23104. / 27. * zeta3 + 20320. / 81. * zeta2 - 24137. / 81. )
            + nf * CA * CF * ( 160. * zeta5 - 128. * zeta3 * zeta2 - 352. / 5. * zeta22
                               + 3712. / 9. * zeta3 + 440. / 3. * zeta2 - 34066. / 81. )
            + nf * CF2 * ( - 320. * zeta5 + 592. / 3. * zeta3 + 572. / 9. )
            + nf * dFF4nc / CF * ( - 1280. / 3. * zeta5 - 256. / 3. * zeta3 + 256. * zeta2 )
            + dFA4nc / CF * ( - 384. * zeta32 - 7936. / 35. * zeta23
                              + 3520. / 3. * zeta5 + 128. / 3. * zeta3 - 128. * zeta2 )
            + CA3 * ( - 16. * zeta32 - 20032. / 105. * zeta23 - 3608. / 9. * zeta5
                      - 352. / 3. * zeta3 * zeta2 + 3608. / 5. * zeta22
                      + 20944. / 27. * zeta3 - 88400. / 81. * zeta2 + 84278. / 81. ) );
  }

  //_________________________________________________________________________
  double gammaK3gmq(int const& nf)
  {
    return
      2 * ( nf * 5. / 24. * ( - 1280. / 3. * zeta5 - 256. / 3. * zeta3 + 256. * zeta2 )
            + 15. / 4. * ( - 384. * pow(zeta3, 2) - 7936. / 35. * pow(zeta2, 3)
                           + 3520. / 3. * zeta5 + 128. / 3. * zeta3 - 128. * zeta2 ) );
  }
}
