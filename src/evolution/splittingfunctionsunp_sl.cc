//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionsunp_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  P0ns::P0ns():
    Expression()
  {
  }
  double P0ns::Regular(double const& x) const
  {
    return - 2 * CF * ( 1 + x );
  }
  double P0ns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double P0ns::Local(double const& x) const
  {
    return 4 * CF * log(1 - x) + 3 * CF;
  }

  //_________________________________________________________________________________
  P0qg::P0qg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0qg::Regular(double const& x) const
  {
    return 4 * _nf * TR * ( 1 - 2 * x + 2 * x * x );
  }

  //_________________________________________________________________________________
  P0gq::P0gq():
    Expression()
  {
  }
  double P0gq::Regular(double const& x) const
  {
    return 4 * CF * ( - 1 + 0.5 * x + 1 / x );
  }

  //_________________________________________________________________________________
  P0gg::P0gg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0gg::Regular(double const& x) const
  {
    return 4 * CA * ( - 2 + x - x * x + 1 / x );
  }
  double P0gg::Singular(double const& x) const
  {
    return 4 * CA / ( 1 - x );
  }
  double P0gg::Local(double const& x) const
  {
    return 4 * CA * log(1 - x) - 2 / 3. * _nf + 11 / 3. * CA;
  }

  //_________________________________________________________________________________
  P1nsp::P1nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = - 40 / 9. * CF * _nf + ( 268 / 9. - 8 * zeta2 ) * CA * CF;
  }
  double P1nsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double gqq1  =
      + 2 * CF * _nf * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 )
      + 4 * CA * CF * ( ( 67 / 18. + 11 * lnx / 6 + lnx2 / 2 - Pi2 / 6 ) * pqq
                        + 20 * ( 1 - x ) / 3 + lnx * ( 1 + x ) )
      + 4 * CF * CF * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x )
                        - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
      + 4 * CF * ( CF - CA / 2 ) * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    return gqq1 - gqq1l;
  }
  double P1nsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1nsp::Local(double const& x) const
  {
    const double p1delta =
      - 1 / 3. * CF * _nf + 3 / 2. * CF * CF + 17 / 6. * CA * CF + 24 * zeta3 * CF * CF - 12 * zeta3 * CA * CF
      - 8 / 3. * zeta2 * CF * _nf - 12 * zeta2 * CF * CF + 44 / 3. * zeta2 * CA * CF;
    return log(1 - x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P1nsm::P1nsm(int const& nf):
    P1nsp(nf)
  {
  }
  double P1nsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double gqq1  =
      + 2 * CF * _nf * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 )
      + 4 * CA * CF * ( ( 67 / 18. + 11 * lnx / 6 + lnx2 / 2 - Pi2 / 6 ) * pqq
                        + 20 * ( 1 - x ) / 3 + lnx * ( 1 + x ) )
      + 4 * CF * CF * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x )
                        - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
      - 4 * CF * ( CF - CA / 2 ) * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    return gqq1 - gqq1l;
  }

  //_________________________________________________________________________________
  P1ps::P1ps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1ps::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return
      _nf * CF * ( - 8 + 24 * x - 224 / 9. * x * x + 80 / 9. / x
                   + ( 4 + 20 * x ) * lnx + 32 / 3. * x * x * lnx
                   - ( 4 + 4 * x ) * lnx2 );
  }

  //_________________________________________________________________________________
  P1qg::P1qg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1qg::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqg   = x * x + ( 1 - x ) * ( 1 - x );
    const double pqgmx = x * x + ( 1 + x ) * ( 1 + x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 2 * CF * _nf * ( 4  + 4 * ln1mx + ( 10 - 4 * ( ln1mx - lnx ) + 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 2 * Pi2 / 3 ) * pqg
                         - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) - 9 * x )
      + 2 * CA * _nf * ( 182 / 9. - 4 * ln1mx
                         + ( - 218 / 9. + 4 * ln1mx - 2 * ln1mx * ln1mx + 44 * lnx / 3 - lnx2 + Pi2 / 3 ) * pqg
                         + 2 * pqgmx * S2x + 40 / ( 9 * x ) + 14 * x / 9 - lnx2 * ( 2 + 8 * x )
                         + lnx * ( - 38 / 3. + 136 * x / 3 ) );
  }

  //_________________________________________________________________________________
  P1gq::P1gq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1gq::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pgq   = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
    const double pgqmx = - ( 1 + ( 1 + x ) * ( 1 + x ) ) / x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 2 * CF * _nf * ( - ( 20 / 9. + 4 * ln1mx / 3 ) * pgq - 4 * x / 3 )
      + 4 * CF * CF  * ( - 2.5 - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq - lnx2 * ( 1 - x / 2 ) - 7 * x / 2
                         - 2 * ln1mx * x + lnx * ( 2 + 7 * x / 2 ) )
      + 4 * CA * CF  * ( 28 / 9. + pgq * ( 0.5 + 11 * ln1mx / 3 + ln1mx * ln1mx - 2 * ln1mx * lnx + lnx2 / 2 - Pi2 / 6 ) + pgqmx * S2x
                         + 65 * x / 18 + 2 * ln1mx * x + 44 * x * x / 9 + lnx2 * ( 4 + x ) - lnx * ( 12 + 5 * x + 8 * x * x / 3 ) );
  }

  //_________________________________________________________________________________
  P1gg::P1gg(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2g = - 40 / 9. * CA * _nf + ( 268 / 9. - 8 * zeta2 ) * CA * CA;
  }
  double P1gg::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pgg   = ( 1 / ( 1 - x ) +  1 / x - 2 + x * ( 1 - x ) );
    const double pggmx = ( 1 / ( 1 + x ) -  1 / x - 2 - x * ( 1 + x ) );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double ggg1  =
      + 2 * CF * _nf * ( - 16 + 4 / ( 3 * x ) + 8 * x + ( 20 * x * x ) / 3 - lnx2 * ( 2 + 2 * x ) - lnx * ( 6 + 10 * x ) )
      + 2 * CA * _nf * ( 2 - 20 * pgg / 9 - 2 * x - 4 * lnx * ( 1 + x ) / 3 + 26 * ( - 1 / x + x * x ) / 9 )
      + 4 * CA *  CA * ( pgg * ( 67 / 9. - 4 * ln1mx * lnx + lnx2 - Pi2 / 3 ) + 2 * pggmx * S2x
                         + 27 * ( 1 - x ) / 2 + 4 * lnx2 * ( 1 + x ) + 67 * ( - 1 / x + x * x ) / 9
                         - lnx * ( 25 / 3. - 11 * x / 3 + 44 * x * x / 3 ) );
    const double ggg1l = _a2g / ( 1 - x );
    return ggg1 - ggg1l;
  }
  double P1gg::Singular(double const& x) const
  {
    return _a2g / ( 1 - x );
  }
  double P1gg::Local(double const& x) const
  {
    const double p1delta = ( - 2 * CF - 8 / 3. * CA ) * _nf + ( 32 / 3. + 12 * zeta3 ) * CA * CA;
    return log(1 - x) * _a2g + p1delta;
  }

  //_________________________________________________________________________________
  P2nsp::P2nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2nsp::Regular(double const& x) const
  {
    const double x_2  = x * x;
    const double x_3  = x_2 * x;
    const double dl   = log(x);
    const double dl_2 = dl * dl;
    const double dl_3 = dl_2 * dl;
    const double dl_4 = dl_3 * dl;
    const double dl1  = log(1 - x);
    const double d81  = 1. / 81.;
    return
      1641.1 - 3135. * x + 243.6 * x_2 - 522.1 * x_3
      + 128. * d81 * dl_4 + 2400. * d81 * dl_3
      + 294.9 * dl_2 + 1258. * dl
      + 714.1 * dl1 + dl * dl1 * ( 563.9 + 256.8 * dl )
      + _nf * ( -197.0 + 381.1 * x + 72.94 * x_2 + 44.79 * x_3
                - 192. * d81 * dl_3  - 2608. * d81 * dl_2 - 152.6 * dl
                - 5120. * d81 * dl1 - 56.66 * dl * dl1 - 1.497 * x * dl_3 )
      + _nf * _nf * ( 32. * x * dl / ( 1 - x ) * ( 3. * dl + 10. ) + 64.
                      + ( 48. * dl_2 + 352. * dl + 384. ) * ( 1 - x ) ) * d81;
  }
  double P2nsp::Singular(double const& x) const
  {
    return ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1 - x );
  }
  double P2nsp::Local(double const& x) const
  {
    const double dl1 = log(1 - x);
    return
      1174.898 * dl1 + 1295.624 - 0.24
      - _nf * ( 183.187 * dl1 + 173.938 - 0.011 )
      + _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
  }

  //_________________________________________________________________________________
  P2nsm::P2nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2nsm::Regular(double const& x) const
  {
    const double x_2  = x * x;
    const double x_3  = x_2 * x;
    const double dl   = log(x);
    const double dl_2 = dl * dl;
    const double dl_3 = dl_2 * dl;
    const double dl1  = log(1 - x);
    const double d81  = 1. / 81.;
    return
      1860.2 - 3505.* x + 297.0 * x_2 - 433.2 * x_3
      + 116. * d81 * dl_3 * dl + 2880. * d81 * dl_3
      + 399.2 * dl_2 + 1465.2 * dl
      + 714.1 * dl1 + dl * dl1 * ( 684.0 + 251.2 * dl )
      + _nf * ( -216.62 + 406.5 * x + 77.89 * x_2 + 34.76 * x_3
                - 256. * d81 * dl_3  - 3216. * d81 * dl_2 - 172.69 * dl
                - 5120. * d81 * dl1 - 65.43 * dl * dl1 - 1.136 * x * dl_3 )
      + _nf * _nf * ( 32.* x * dl / ( 1 - x ) * ( 3. * dl + 10. ) + 64.
                      + ( 48.* dl_2 + 352.* dl + 384. ) * ( 1.-x ) ) * d81;
  }
  double P2nsm::Singular(double const& x) const
  {
    return ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1 - x );
  }
  double P2nsm::Local(double const& x) const
  {
    const double dl1 = log(1 - x);
    return
      1174.898 * dl1 + 1295.624 - 0.154
      - _nf * ( 183.187 * dl1 + 173.938  - 0.005 )
      + _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
  }

  //_________________________________________________________________________________
  P2nss::P2nss(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2nss::Regular(double const& x) const
  {
    const double x_2  = x * x;
    const double d27  = 1. / 27.;
    const double dl   = log(x);
    const double dl_2 = dl * dl;
    const double dl_3 = dl_2 * dl;
    const double dl_4 = dl_3 * dl;
    const double x1   = 1 - x;
    const double dl1  = log(x1);
    return
      _nf * ( x1 * ( 151.49 + 44.51 * x - 43.12 * x_2 + 4.820 * x_2 * x )
              + 40. * d27 * dl_4 - 80. * d27 * dl_3 + 6.892 * dl_2
              + 178.04 * dl + dl * dl1 * ( - 173.1 + 46.18 * dl )
              + x1 * dl1 * ( - 163.9 / x - 7.208 * x ) );
  }

  //_________________________________________________________________________________
  P2ps::P2ps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2ps::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1 - x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double  p2ps1 =
      - 3584. / ( 27. * x ) * dl - 506.0 / x + 160. / 27. * dl_4
      - 400. / 9. * dl_3 + 131.4 * dl_2 - 661.6 * dl
      - 5.926 * dl1_3 - 9.751 * dl1_2 - 72.11 * dl1
      + 177.4 + 392.9 * x - 101.4 * x_2 - 57.04 * dl * dl1;
    const double p2ps2  =
      256. / ( 81. * x ) + 32. / 27. * dl_3 + 17.89 * dl_2
      + 61.75 * dl + 1.778 * dl1_2 + 5.944 * dl1 + 100.1
      - 125.2 * x + 49.26 * x_2 - 12.59 * x_3
      - 1.889 * dl * dl1;
    return ( 1 - x ) * _nf * ( p2ps1 + _nf * p2ps2 );
  }

  //_________________________________________________________________________________
  P2qg::P2qg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2qg::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1 - x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double dl1_4 = dl1_3 * dl1;
    const double p2qg1 =
      - 896. / ( 3. * x ) * dl - 1268.3 / x + 536./27. * dl_4
      - 44. / 3. * dl_3 + 881.5 * dl_2 + 424.9 * dl
      + 100. / 27. * dl1_4 - 70. / 9. * dl1_3
      - 120.5 * dl1_2 + 104.42 * dl1
      + 2522. - 3316. * x + 2126. * x_2
      + dl * dl1 * ( 1823. - 25.22 * dl ) - 252.5 * x * dl_3;
    const double p2qg2 =
      1112. / ( 243. * x ) - 16. / 9. * dl_4
      - 376. / 27. * dl_3 - 90.8 * dl_2 - 254.0 * dl
      + 20./27. * dl1_3 + 200. / 27. * dl1_2 - 5.496 * dl1
      - 252.0  + 158.0 * x + 145.4 * x_2 - 139.28 * x_3
      - dl * dl1 * ( 53.09  + 80.616 * dl ) - 98.07 * x * dl_2
      + 11.70 * x * dl_3;
    return _nf * ( p2qg1 + _nf * p2qg2 );
  }

  //_________________________________________________________________________________
  P2gq::P2gq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2gq::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1 - x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double dl1_4 = dl1_3 * dl1;
    const double p2gq0 =
      1189.3 * dl / x + 6163.1 / x - 4288. / 81. * dl_4
      + 1568. / 9. * dl_3 - 1794. * dl_2 + 4033. * dl
      + 400. / 81. * dl1_4 + 2200. / 27. * dl1_3
      + 606.3 * dl1_2 + 2193. * dl1
      - 4307. + 489.3 * x + 1452.* x_2 + 146.0 * x_3
      - 447.3 * dl_2 * dl1 - 972.9 * x * dl_2;
    const double p2gq1 =
      71.082 * dl / x  - 46.41 / x + 128. / 27. * dl_4
      + 704/81. * dl_3 + 20.39 * dl_2 + 174.8 * dl
      - 400./81. * dl1_3 - 68.069 * dl1_2 - 296.7 * dl1
      - 183.8 + 33.35 * x - 277.9 * x * x + 108.6 * x * dl_2
      - 49.68 * dl * dl1;
    const double p2gq2 =
      ( 64. * ( - 1. / x + 1. + 2. * x )
        + 320. * dl1 * ( 1. / x - 1. + 0.8 * x )
        + 96. * dl1_2 * ( 1. / x - 1. + 0.5 * x ) ) / 27.;
    return ( p2gq0 + _nf * ( p2gq1 + _nf * p2gq2 ) );
  }

  //_________________________________________________________________________________
  P2gg::P2gg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2gg::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1 - x);
    const double p2gga0 =
      2675.8 * dl / x + 14214. / x - 144. * dl_4 + 72. * dl_3
      - 7471. * dl_2 + 274.4 * dl + 3589. * dl1 - 20852.
      + 3968. * x - 3363. * x_2 + 4848. * x_3
      + dl * dl1 * ( 7305. + 8757. * dl );
    const double p2gga1 =
      157.27 * dl / x + 182.96 / x + 512./27. * dl_4
      + 832. / 9. * dl_3 + 491.3 * dl_2 + 1541. * dl
      - 320.0 * dl1 - 350.2 + 755.7 * x - 713.8 * x_2
      + 559.3 * x_3 + dl * dl1 * ( 26.15 - 808.7 * dl );
    const double p2gga2 =
      - 680. / ( 243. * x ) - 32. / 27. * dl_3 + 9.680 * dl_2
      - 3.422 * dl - 13.878 + 153.4 * x - 187.7 * x_2
      + 52.75 * x_3 - dl * dl1 * ( 115.6 - 85.25 * x + 63.23 * dl);
    return p2gga0 + _nf * ( p2gga1 + _nf * p2gga2 );
  }
  double P2gg::Singular(double const& x) const
  {
    return ( 2643.521 - _nf * 412.172 - _nf * _nf * 16. / 9. ) / ( 1 - x );
  }
  double P2gg::Local(double const& x) const
  {
    const double dl1 = log(1 - x);
    return
      2643.521 * dl1 + 4425.448 + 0.446
      - _nf * ( 412.172 * dl1 +  528.720 + 0.003 )
      + _nf * _nf * ( - 16. / 9. * dl1 + 6.4630 );
  }

  // //_________________________________________________________________________________
  // P3nsp::P3nsp(int const& nf, int const& imod, double const& rho):
  //   Expression(),
  //   _nf(nf),
  //   _imod(imod),
  //   _rho(rho)
  // {
  //   const int nf2 = _nf * _nf;
  //   const int nf3 = _nf * nf2;

  //   // Moments for the known exact small-x and large-x contributions (Vogt)
  //   std::vector<double> N(8, 0.);
  //   N[0] = -0.00021211
  //          - _nf * ( 0.1654801313006386 + 0.20740623526104135 - 0.1502823524 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( -0.004654314166031352 -0.007842937752425193 + 4.7723213922e-03 + 7.774039377e-03 )
  //          - nf3 * ( -6.62104090137594e-05 - 0.0001312258806063389 + 7.653104787e-05 + 1.209052801e-04 )
  //          + _rho * 0.375;
  //   N[1] = -0.00021141
  //          - _nf * ( 0.11120681844515805 + 0.3802447646452425 - 0.2688473229 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( 9.571697043e-3 + 7.774039377e-03 -0.0029670140315604917 - 0.014378719212779521 )
  //          - nf3 * ( -3.609651149299798e-05 -0.00024058078111162135 + 1.55772043e-04 + 1.209052801e-04 )
  //          + _rho * 0.0234375;
  //   N[2] = -0.00021069
  //          - _nf * ( 0.08533032993183774 + 0.4735775705127111 - 0.3363034306 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( 0.01237066021 + 7.774039377e-03 - 0.002236654071582081 -0.01790804120137086 )
  //          - nf3 * ( -2.4727034394311423e-05 -0.0002996324273844739 + 2.034542201e-04 + 1.209052801e-04 )
  //          + _rho * 0.004629629629528962;
  //   N[3] = -0.00020986
  //          - _nf * ( 0.06978655862237221 + 0.5377747385697003 - 0.3849568093 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( -0.001806806533621782 -0.020335617172359613 + 0.01436839622 + 7.774039377e-03 )
  //          - nf3 * ( -1.875753442493042e-05 -0.00034024996185786445 + 2.381022604e-04 + 1.209052801e-04 )
  //          + _rho * 0.0014648437499999065;
  //   N[4] = -0.00020809
  //          - _nf * ( 0.059311724466073266 + 0.5867456552285572 - 0.4234530704 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( -0.0015200868207649136 -0.02218742191946 + 0.01593349271 + 7.774039377e-03 )
  //          - nf3 * ( -1.5092063174265756e-05 - 0.0003712338503343612 + 2.654206822e-04 + 1.209052801e-04 )
  //          + _rho * 0.0005999999999999996;
  //   N[5] = -0.0002073
  //          - _nf * ( 0.0517320184770643 + 0.6263413910511195 - 0.4554694245 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( -0.001314205205204153 -0.023684710035832085 + 0.0172249122 + 7.774039377e-03 )
  //          - nf3 * ( -1.2617371403920047e-05 -0.0003962860639046623 + 2.879982076e-04 + 1.209052801e-04 )
  //          + _rho * 0.0002893518518518518;
  //   N[6] = -0.00020562
  //          - _nf * ( 0.045972268757607376 + 0.6595795697788507 - 0.4829482701 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( -0.0011588022146498728 -0.024941591085900227 + 0.0183264033 + 7.774039377e-03 )
  //          - nf3 * ( -1.0836262424080888e-05 -0.0004173158524633704 + 3.072468902e-04 + 1.209052801e-04 )
  //          + _rho * 0.00015618492294877136;
  //   N[7] = -0.00020381
  //          - _nf * ( 0.04143575048971359 + 0.6882213832196611 - 0.5070540331 + 0.01076711056 - 0.2333500675 )
  //          - nf2 * ( -0.0010371493846056854 -0.026024663442187512 + 0.01928783511 + 7.774039377e-03 )
  //          - nf3 * ( -9.493881167514427e-06 -0.00043543752168996016 + 3.240261808e-04 + 1.209052801e-04 )
  //          + _rho * 9.155273437499997e-05;

  //   // Matrix
  //   const std::vector<std::vector<double>>  inv_A
  //   {
  //     {9.67452285e-02, 7.12205880e03, 2.55568773e03, 4.41612871e02, -3.63251998e02, 8.54599881e03,  -8.03082218e03, 5.66780493e01},
  //     {-1.10377154e01, -7.11199538e05, -2.58657431e05, -4.59602884e04, 2.63108547e04, -8.24798082e05, 7.81924565e05, -3.25280541e03},
  //     {1.99489100e02, 1.15489637e07, 4.23773320e06, 7.69318756e05, -3.54975874e05, 1.31604354e07, -1.25178270e07, 4.09448382e04},
  //     {-1.28368247e03, -6.79475656e07, -2.50856870e07, -4.63268153e06, 1.85137884e06, -7.66032630e07, 7.29666144e07, -2.06388099e05},
  //     {3.82370653e03, 1.87306734e08, 6.94570801e07, 1.30090317e07, -4.68326682e06, 2.09632375e08, -1.99805999e08, 5.11458122e05},
  //     {-5.75484319e03, -2.63212244e08, -9.79232841e07, -1.85597818e07, 6.16625773e06, -2.93018936e08, 2.79351504e08, -6.64139088e05},
  //     {4.25513991e03, 1.82949771e08, 6.82322632e07, 1.30646558e07, -4.07104678e06, 2.02835925e08, -1.93381943e08, 4.34087442e05},
  //     {-1.22905233e03, -4.99408525e07, -1.86617424e07, -3.60495568e06, 1.06569809e06, -5.51895151e07, 5.26130529e07, -1.12766384e05}
  //   };

  //   // Matrix multiplication of inv_A and N
  //   _C.resize(N.size(), 0.);
  //   for (int i = 0; i < (int) N.size(); i++)
  //     for (int j = 0; j < (int) N.size(); j++)
  //       _C[i] += inv_A[j][i] * N[j];
  // }
  // double P3nsp::Regular(double const& y) const
  // {
  //   const double y2   = y * y;
  //   const double y3   = y2 * y;
  //   const double omy  = 1 - y;
  //   const double dm   = 1 / omy;
  //   const double dl   = log(y);
  //   const double dl2  = dl * dl;
  //   const double dl3  = dl2 * dl;
  //   const double dl4  = dl3 * dl;
  //   const double dl5  = dl4 * dl;
  //   const double dl6  = dl5 * dl;
  //   const double dlm  = log(omy);
  //   const double dlm2 = dlm * dlm;
  //   const double dlm3 = dlm2 * dlm;

  //   // Leading large-n_c, nf^0 and nf^1, parametrized.
  //   const double p3nsa0 =
  //     2.5e4 * ( omy * ( 3.5254 + 8.6935 * y - 1.5051 * y2 + 1.8300 * y3 )
  //               + 11.883 * y * dl - 0.09066 * y * dl2 + 11.410 * omy * dlm + 13.376 * dl * dlm )
  //     + 5.167133e4 * dl + 1.712095e4 * dl2 + 2.863226e3 * dl3 + 2.978255e2 * dl4
  //     + 1.6e1 * dl5 + 5e-1 * dl6 - 2.973385e4 + 1.906980e4 * dlm;
  //   const double p3nsa1 =
  //     2.5e4* ( omy * ( - 0.74077 + 1.4860 * y - 0.23631 * y2 + 0.31584 * y3 )
  //              + 2.5251 * omy * dlm + 2.5203 * dl * dlm + 2.2242 * y * dl
  //              - 0.02460 * y * dl2 + 0.00310 * y * dl3 )
  //     - 9.239374e3 * dl - 2.917312e3 * dl2 - 4.305308e2 * dl3 - 3.6e1 * dl4
  //     - 4. / 3. * dl5 + 8.115605e3 - 3.079761e3 * dlm;

  //   // Nonleading large-n_c, nf^0 and nf^1: two approximations
  //   const double p3npa01 =
  //     3948.16 * omy - 2464.61 * ( 2 * y - y2 ) * omy - 1839.44 * dl2 - 402.156 * dl3
  //     - 1777.27 * dlm2 * omy - 204.183 * dlm3 * omy + 507.152 - 5.587553e+1 * dl4
  //     - 2.831276 * dl5 - 1.488340e-1 * dl6 - 2.601749e+3 - 2.118867e+3 * dlm;
  //   const double p3npa02 =
  //     ( 8698.39 - 10490.47 * y ) * y * omy + 1389.73 * dl + 189.576 * dl2
  //     - 173.936 * dlm2 * omy + 223.078 * dlm3 * omy + 505.209 - 5.587553e+1 * dl4
  //     - 2.831276 * dl5 - 1.488340e-1 * dl6 - 2.601749e+3 - 2.118867e+3 * dlm;

  //   const double p3npa11 =
  //     ( - 1116.34 + 1071.24 * y ) * y * omy - 59.3041 * dl2 - 8.4620 * dl3
  //     - 143.813 * dlm * omy - 18.8803 * dlm3 * omy - 7.33927 + 4.658436 * dl4
  //     + 2.798354e-1 * dl5 + 3.121643e+2 + 3.379310e+2 * dlm;
  //   const double p3npa12 =
  //     ( - 690.151 - 656.386 * y2 ) * omy + 133.702 * dl2 + 34.0569 * dl3
  //     - 745.573 * dlm * omy + 8.61438 * dlm3 * omy - 7.53662 + 4.658437 * dl4
  //     + 2.798354e-1 * dl5 + 3.121643e+2 + 3.379310e+2 * dlm;

  //   // nf^2 (parametrized) and nf^3 (exact)
  //   const double p3nspa2 =
  //     2.5e2 *  ( omy * ( 3.0008 + 0.8619 * y - 0.12411 * y2 + 0.31595 * y3 )
  //                - 0.37529 * y * dl - 0.21684 * y * dl2 - 0.02295 * y * dl3
  //                + 0.03394 * omy * dlm + 0.40431 * dl * dlm )
  //     + 3.930056e+2 * dl + 1.125705e+2 * dl2 + 1.652675e+1 * dl3
  //     + 7.901235e-1 * dl4 - 3.760092e+2 + 2.668861e+1 * dlm;
  //   const double p3nsa3  =
  //     - 2.426296e0 - 8.460488e-1 * y + ( 5.267490e-1 * dm - 3.687243e0 + 3.160494e0 * y ) * dl
  //     - ( 1.316872e0 * ( dm + 1e-1 ) - 1.448560e0 * y ) * dl2
  //     - ( 2.633745e-1 * dm - 1.31687e-1 * ( 1 + y ) ) * dl3;

  //   // Assembly
  //   const double p3nspai = p3nsa0 + _nf * ( p3nsa1 + _nf * ( p3nspa2 + _nf * p3nsa3 ) )
  //                          + ( _C[1] * omy * dlm + _C[2] * omy * dlm2 + _C[3] * omy * dlm3
  //                              + _C[4] + _C[5] * y + _C[6] * y2 + _C[7] * dl2 + _rho * dl3 ) * pow(FourPi, 4);
  //   if (_imod == 1)
  //     return p3nspai + p3npa01 + _nf * p3npa11;
  //   else if (_imod == 2)
  //     return p3nspai + p3npa02 + _nf * p3npa12;
  //   else
  //     return p3nspai + 0.5 * ( ( p3npa01 + p3npa02 ) + _nf * ( p3npa11 + p3npa12 ) );
  // }
  // double P3nsp::Singular(double const& y) const
  // {
  //   const double d1 = 1 / ( 1 - y );
  //   const double a4qi =
  //     2.120902e+4
  //     - 5.179372e+3 * _nf
  //     + 1.955772e+2 * _nf * _nf
  //     + 3.272344 * _nf * _nf * _nf
  //     + _C[0] * pow(FourPi, 4);
  //   const double a4ap1 = - 507.152 + 7.33927 * _nf;
  //   const double a4ap2 = - 505.209 + 7.53662 * _nf;

  //   if (_imod == 1)
  //     return ( a4qi + a4ap1 ) * d1;
  //   else if (_imod == 2)
  //     return ( a4qi + a4ap2 ) * d1;
  //   else
  //     return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * d1;
  // }
  // double P3nsp::Local(double const& y) const
  // {
  //   const double dl1 = log(1 - y);
  //   const double a4qi  =
  //     2.120902e+4
  //     - 5.179372e+3 * _nf
  //     + 1.955772e+2 * _nf * _nf
  //     + 3.272344 * _nf * _nf * _nf
  //     + _C[0] * pow(FourPi, 4);
  //   const double a4ap1 = - 507.152 + 7.33927 * _nf;
  //   const double a4ap2 = - 505.209 + 7.53662 * _nf;

  //   const double b4qi =
  //     2.579609e+4 + 0.08
  //     - ( 5.818637e+3 + 0.97 ) * _nf
  //     + ( 1.938554e+2 + 0.0037)* _nf * _nf
  //     + 3.014982 * _nf * _nf * _nf;
  //   const double b4ap1 = - 2405.03 + 267.965 * _nf;
  //   const double b4ap2 = - 2394.47 + 269.028 * _nf;

  //   if (_imod == 1)
  //     return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap1;
  //   else if (_imod == 2)
  //     return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap2;
  //   else
  //     return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * dl1
  //            + b4qi + 0.5 * ( b4ap1 + b4ap2 );
  // }

  //_________________________________________________________________________________
  P3nsp::P3nsp(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double P3nsp::Regular(double const& y) const
  {
    const double y2   = y * y;
    const double y3   = y2 * y;
    const double omy  = 1 - y;
    const double dm   = 1 / omy;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl2 * dl;
    const double dl4  = dl3 * dl;
    const double dl5  = dl4 * dl;
    const double dl6  = dl5 * dl;
    const double dlm  = log(omy);
    const double dlm2 = dlm * dlm;
    const double dlm3 = dlm2 * dlm;

    // Leading large-n_c, nf^0 and nf^1, parametrized
    const double p3nsa0  =
      2.5e+4 * ( omy * ( 3.5254 + 8.6935 * y - 1.5051 * y2 + 1.8300 * y3 )
                 + 11.883 * y * dl - 0.09066 * y * dl2 + 11.410 * omy * dlm + 13.376 * dl * dlm )
      + 5.167133e+4 * dl + 1.712095e+4 * dl2 + 2.863226e+3 * dl3 + 2.978255e+2 * dl4
      + 1.6e+1 * dl5 + 5.e-1 * dl6 - 2.973385e+4 + 1.906980e+4 * dlm;
    const double p3nsa1  =
      2.5e+4 * ( omy * ( - 0.74077 + 1.4860 * y - 0.23631 * y2 + 0.31584 * y3 )
                 + 2.5251 * omy * dlm + 2.5203 * dl * dlm + 2.2242 * y * dl
                 - 0.02460 * y * dl2 + 0.00310 * y * dl3 )
      - 9.239374e+3 * dl - 2.917312e+3 * dl2 - 4.305308e+2 *dl3 - 3.6e+1 * dl4
      - 4. / 3. * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

    // Nonleading large-n_c, nf^0 and nf^1: two approximations
    const double p3npa01 =
      3948.16 * omy - 2464.61 * ( 2 * y - y2 ) * omy - 1839.44 * dl2 - 402.156 * dl3
      - 1777.27 * dlm2 * omy - 204.183 * dlm3 * omy + 507.152 - 5.587553e1 * dl4 - 2.831276e0 * dl5
      - 1.488340e-1 * dl6 - 2.601749e3 - 2.118867e3 * dlm;
    const double p3npa02 =
      ( 8698.39 - 10490.47 * y ) * y * omy + 1389.73 * dl + 189.576 * dl2
      - 173.936 * dlm2 * omy + 223.078 * dlm3 * omy + 505.209 - 5.587553e1 * dl4 - 2.831276e0 * dl5
      - 1.488340e-1 * dl6 - 2.601749e3 - 2.118867e3 * dlm;

    const double p3npa11 =
      ( - 1116.34 + 1071.24 * y ) * y * omy - 59.3041 * dl2 - 8.4620 * dl3
      - 143.813 * dlm * omy - 18.8803 * dlm3 * omy - 7.33927 + 4.658436e0*dl4 + 2.798354e-1 * dl5
      + 3.121643e2 + 3.379310e2 * dlm;
    const double p3npa12 =
      ( - 690.151 - 656.386 * y2 ) * omy + 133.702 * dl2 + 34.0569 * dl3
      - 745.573 * dlm * omy + 8.61438 * dlm3 * omy - 7.53662 + 4.658437e0 * dl4 + 2.798354e-1 * dl5
      + 3.121643e2 + 3.379310e2 * dlm;

    // nf^2 (parametrized) and nf^3 (exact)
    const double p3nspa2 =
      2.5e+2 * ( omy * ( 3.0008 + 0.8619 * y - 0.12411 * y2 + 0.31595* y3 )
                 - 0.37529 * y * dl - 0.21684 * y * dl2 - 0.02295 * y * dl3
                 + 0.03394 * omy * dlm + 0.40431 * dl * dlm )
      + 3.930056e+2 * dl + 1.125705e+2 * dl2 + 1.652675e+1 * dl3
      + 7.901235e-1 * dl4 - 3.760092e+2 + 2.668861e+1 * dlm;
    const double p3nsa3  =
      - 2.426296 - 8.460488e-1 * y + ( 5.267490e-1 * dm - 3.687243 + 3.160494 * y ) * dl
      - ( 1.316872 * ( dm + 1e-1) - 1.448560 * y ) * dl2
      - ( 2.633745e-1 * dm - 1.31687e-1 * ( 1 + y ) ) * dl3;

    // Assembly
    const double p3nspai = p3nsa0 + _nf * p3nsa1 + _nf * _nf * p3nspa2 + _nf * _nf * _nf * p3nsa3;
    if (_imod == 1)
      return p3nspai + p3npa01 + _nf * p3npa11;
    else if (_imod == 2)
      return p3nspai + p3npa02 + _nf * p3npa12;
    else
      return p3nspai + 0.5 * ( ( p3npa01 + p3npa02 ) + _nf * ( p3npa11 + p3npa12 ) );
  }
  double P3nsp::Singular(double const& y) const
  {
    const double d1 = 1 / ( 1 - y );

    const double a4qi  =
      2.120902e+4
      - 5.179372e+3 * _nf
      + 1.955772e+2 * _nf * _nf
      + 3.272344e+0 * _nf * _nf * _nf;
    const double a4ap1 = - 507.152 + 7.33927 * _nf;
    const double a4ap2 = - 505.209 + 7.53662 * _nf;

    if (_imod == 1)
      return ( a4qi + a4ap1 ) * d1;
    else if (_imod == 2)
      return ( a4qi + a4ap2 ) * d1;
    else
      return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * d1;
  }
  double P3nsp::Local(double const& y) const
  {
    const double dl1 = log(1 - y);

    const double a4qi  =
      2.120902e+4
      - 5.179372e+3 * _nf
      + 1.955772e+2 * _nf * _nf
      + 3.272344e+0 * _nf * _nf * _nf;
    const double a4ap1 = - 507.152 + 7.33927 * _nf;
    const double a4ap2 = - 505.209 + 7.53662 * _nf;

    const double b4qi =
      2.579609e+4 + 0.08
      - ( 5.818637e+3 + 0.97 ) * _nf
      + ( 1.938554e+2 + 0.0037 ) * _nf * _nf
      +   3.014982e+0 * _nf * _nf * _nf;
    const double b4ap1 = - 2405.03 + 267.965 * _nf;
    const double b4ap2 = - 2394.47 + 269.028 * _nf;

    if (_imod == 1)
      return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap1;
    else if (_imod == 2)
      return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap2;
    else
      return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * dl1 + b4qi + 0.5 * ( b4ap1 + b4ap2 );
  }

  //_________________________________________________________________________________
  P3nsm::P3nsm(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double P3nsm::Regular(double const& y) const
  {
    const double y2   = y * y;
    const double y3   = y2 * y;
    const double omy  = 1 - y;
    const double dm   = 1 / omy;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl2 * dl;
    const double dl4  = dl3 * dl;
    const double dl5  = dl4 * dl;
    const double dl6  = dl5 * dl;
    const double dlm  = log(omy);
    const double dlm2 = dlm * dlm;
    const double dlm3 = dlm2 * dlm;

    // Leading large-n_c, nf^0 and nf^1, parametrized
    const double p3nsa0  =
      2.5e+4 * ( omy * ( 3.5254 + 8.6935 * y - 1.5051 * y2 + 1.8300 * y3 )
                 + 11.883 * y * dl - 0.09066 * y * dl2 + 11.410 * omy * dlm + 13.376 * dl * dlm )
      + 5.167133e+4 * dl + 1.712095e+4 * dl2 + 2.863226e+3 * dl3 + 2.978255e+2 * dl4
      + 1.6e+1 * dl5 + 5.e-1 * dl6 - 2.973385e+4 + 1.906980e+4 * dlm;
    const double p3nsa1  =
      2.5e+4 * ( omy * ( - 0.74077 + 1.4860 * y - 0.23631 * y2 + 0.31584 * y3 )
                 + 2.5251 * omy * dlm + 2.5203 * dl * dlm + 2.2242 * y * dl
                 - 0.02460 * y * dl2 + 0.00310 * y * dl3 )
      - 9.239374e+3 * dl - 2.917312e+3 * dl2 - 4.305308e+2 *dl3 - 3.6e+1 * dl4
      - 4. / 3. * dl5 + 8.115605e+3 - 3.079761e+3 * dlm;

    // Nonleading large-n_c, nf^0 and nf^1: two approximations
    const double p3nma01 =
      ( 5992.88 * ( 1 + 2 * y ) + 31321.44 * y2 ) * omy + 511.228 - 1618.07 * dl + 2.25480 * dl3
      + 31897.82 * dlm * omy + 4653.76 * dlm2 * omy + 4.964335e-1 * ( dl6 + 6 * dl5 )
      - 2.601749e+3 - 2.118867e+3 * dlm;
    const double p3nma02 =
      ( 4043.59 - 15386.6 * y ) * y * omy + 502.481 + 1532.96 * dl2 + 31.6023 * dl3
      - 3997.39 * dlm * omy + 511.567 * dlm3 * omy + 4.964335e-1 * ( dl6 + 18 * dl5 )
      - 2.601749e+3 - 2.118867e+3 * dlm;

    const double p3nma11 =
      ( 114.457 * ( 1 + 2 * y ) + 2570.73 * y2 ) * omy - 7.08645 - 127.012 * dl2 + 2.69618 * dl4
      + 1856.63 * dlm * omy + 440.17 * dlm2 * omy + 3.121643e+2 + 3.379310e+2 * dlm;
    const double p3nma12 =
      ( - 335.995 * ( 2 + y ) - 1605.91 * y2 ) * omy - 7.82077 - 9.76627 * dl2 + 0.14218 * dl5
      - 1360.04 * dlm * omy + 38.7337 * dlm3 * omy + 3.121643e+2 + 3.379310e+2 * dlm;

    // nf^2 (parametrized) and nf^3 (exact)
    const double p3nsma2 =
      2.5e+2 * ( omy * ( 3.2206 + 1.7507 * y + 0.13281 * y2 + 0.45969 * y3 )
                 + 1.5641 * y * dl - 0.37902 * y * dl2 - 0.03248 * y *dl3
                 + 2.7511 * omy * dlm + 3.2709 * dl * dlm )
      + 4.378810e+2 * dl + 1.282948e+2 * dl2 + 1.959945e+1 * dl3
      + 9.876543e-1 * dl4 - 3.760092e+2 + 2.668861e+1 * dlm;
    const double p3nsa3  =
      - 2.426296 - 8.460488e-1 * y + ( 5.267490e-1 * dm - 3.687243 + 3.160494 * y ) * dl
      - ( 1.316872 * ( dm + 1e-1) - 1.448560 * y ) * dl2
      - ( 2.633744e-1 * dm - 1.31687e-1 * ( 1 + y ) ) * dl3;

    // Assembly
    const double p3nsmai = p3nsa0 + _nf * p3nsa1 + _nf * _nf * p3nsma2 + _nf * _nf * _nf * p3nsa3;
    if (_imod == 1)
      return p3nsmai + p3nma01 + _nf * p3nma11;
    else if (_imod == 2)
      return p3nsmai + p3nma02 + _nf * p3nma12;
    else
      return p3nsmai + 0.5 * ( ( p3nma01 + p3nma02 ) + _nf * ( p3nma11 + p3nma12 ) );
  }
  double P3nsm::Singular(double const& y) const
  {
    const double d1 = 1 / ( 1 - y );

    const double a4qi  =
      2.120902e+4
      - 5.179372e+3 * _nf
      + 1.955772e+2 * _nf * _nf
      + 3.272344e+0 * _nf * _nf * _nf;
    const double a4ap1 = - 511.228 + 7.08645 * _nf;
    const double a4ap2 = - 502.481 + 7.82077 * _nf;

    if (_imod == 1)
      return ( a4qi + a4ap1 ) * d1;
    else if (_imod == 2)
      return ( a4qi + a4ap2 ) * d1;
    else
      return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * d1;
  }
  double P3nsm::Local(double const& y) const
  {
    const double dl1 = log(1 - y);

    const double a4qi  =
      2.120902e+4
      - 5.179372e+3 * _nf
      + 1.955772e+2 * _nf * _nf
      + 3.272344e+0 * _nf * _nf * _nf;
    const double a4ap1 = - 511.228 + 7.08645 * _nf;
    const double a4ap2 = - 502.481 + 7.82077 * _nf;

    const double b4qi =
      2.579609e+4 + 0.08
      - ( 5.818637e+3 + 0.97 ) * _nf
      + ( 1.938554e+2 + 0.0037 ) * _nf * _nf
      +   3.014982e+0 * _nf * _nf * _nf;
    const double b4ap1 = - 2426.05  + 266.674 * _nf - 0.05 * _nf;
    const double b4ap2 = - 2380.255 + 270.518 * _nf - 0.05 * _nf;

    if (_imod == 1)
      return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap1;
    else if (_imod == 2)
      return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap2;
    else
      return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * dl1 + b4qi + 0.5 * ( b4ap1 + b4ap2 );
  }

  //_________________________________________________________________________________
  P3nss::P3nss(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double P3nss::Regular(double const& y) const
  {
    const double y2   = y * y;
    const double omy  = 1 - y;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl2 * dl;
    const double dl4  = dl3 * dl;
    const double dl5  = dl4 * dl;
    const double dl6  = dl5 * dl;
    const double dlm  = log(omy);
    const double dlm2 = dlm * dlm;
    const double dlm3 = dlm2 * dlm;

    // nf^1: two approximations
    const double p3nsa11 =
      omy * y * ( 4989.2 - 1607.73 * y ) + 3687.6 * dl + 3296.6 * dl2 + 1271.11* dl3
      + 533.44 * dl4 + 97.27 * dl5 + 4 * dl6 + 60.40 * omy * dlm2 + 4.685 * omy * dlm3;
    const double p3nsa12 =
      1030.79 * omy * y + 1266.77 * omy * ( 2 - y2 ) + 2987.83 * dl + 273.05 * dl2 - 923.48 * dl3
      - 236.76 * dl4 - 33.886 * dl5 - 4 * dl6 - 254.63 * omy * dlm - 0.28953 * omy * dlm3;

    // nf^2 (parametrized)
    const double p3nssa2 =
      2.5e+2 * ( omy * ( - 4.7656 + 1.6908 * y + 0.1703 * y2 )
                 - 0.41652 * y *dl + 0.90777 * y * dl2 + 0.12478 * y * dl3
                 + 0.17155 * omy * dlm + 0.17191 * dl * dlm )
      - 6.473971e+2 * dl - 6.641219e+1 * dl2 - 5.353347 * dl3 - 5.925926 * dl4
      - 3.950617e-1 * dl5 + 1.970002e+1 * omy * dlm - 3.435474 * omy * dlm2;

    if (_imod == 1)
      return _nf * p3nsa11 + _nf * _nf * p3nssa2;
    else if (_imod == 2)
      return _nf * p3nsa12 + _nf * _nf * p3nssa2;
    else
      return 0.5 *_nf * ( p3nsa11 + p3nsa12 ) + _nf * _nf * p3nssa2;
  }

  // //_________________________________________________________________________________
  // P3ps::P3ps(int const& nf, double const& rho):
  //   Expression(),
  //   _nf(nf),
  //   _rho(rho)
  // {
  //   // Moments for the known exact small-x contribution (Vogt)
  //   const std::vector<double> N{-1.28827235e-01 - _rho * 0.75, -2.42108771e-03 - _rho * 0.04861111111078488,
  //                               -6.81270332e-05 - _rho * 0.012222222222363594, 4.79640375e-05 - _rho * 0.004783163265306247};

  //   // Matrix
  //   const std::vector<std::vector<double>>  inv_A
  //   {
  //     {-4.33423515e+01, 1.04394869e+03, -3.26481918e+03, 2.49215628e+03},
  //     {1.01800794e+02, -2.29059693e+03, 7.08126548e+03, -5.37933607e+03},
  //     {1.28895321e+02, -3.28605015e+03, 1.09724714e+04, -8.75612151e+03},
  //     {-9.56650450, 2.57515782e+02, -9.20800163e+02, 7.84949427e+02}
  //   };

  //   // Matrix multiplication of inv_A and N
  //   _C.resize(N.size(), 0.);
  //   for (int i = 0; i < (int) N.size(); i++)
  //     for (int j = 0; j < (int) N.size(); j++)
  //       _C[i] += inv_A[i][j] * N[j];
  // }
  // double P3ps::Regular(double const& x) const
  // {
  //   const double dl  = log(x);
  //   const double dl2 = dl * dl;
  //   const double dlm = log(1 - x);
  //   return _nf * ( ( 1 - x ) * ( _C[0] / x + _C[1] * dl2 + _C[2] * pow(x, 2) + _C[3] * pow(dlm, 2) - _rho * dl / x )
  //                  + pow(CA, 2) * CF * ( 82. / 81. + 2 * zeta3 ) * dl2 / x / 6 / pow(M_PI, 4) ) * pow(FourPi, 4);
  // }

  //_________________________________________________________________________________
  P3ps::P3ps(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double P3ps::Regular(double const& y) const
  {
    const int nf      = _nf;
    const int nf2     = _nf * _nf;
    const int nf3     = _nf * nf2;
    const double ym   = 1 / y;
    const double y1   = 1 - y;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double dl6  = dl * dl5;
    const double dl1  = log(1 - y);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;

    // Known large-x coefficients
    const double x1l4cff = - 5.6460905e1 * nf + 3.6213992e0 * nf2;
    const double x1l3cff = - 2.4755054e2 * nf + 4.0559671e1 * nf2 - 1.5802469e0 * nf3;
    const double y1l4cff = - 1.3168724e+1 * nf;
    const double y1l3cff = - 1.9911111e+2 * nf + 1.3695473e+1 * nf2;

    // Known small-x coefficients
    const double bfkl1   =   1.7492273e3 * nf;
    const double x0l6cff = - 7.5061728e0 * nf + 7.9012346e-1 * nf2;
    const double x0l5cff =   2.8549794e1 * nf + 3.7925926e0 * nf2;
    const double x0l4cff = - 8.5480010e2 * nf + 7.7366255e1 * nf2 - 1.9753086e-1 * nf3;

    // The resulting part of the function
    const double p3ps01 =
      + bfkl1 * dl2 * ym
      + x0l6cff * dl6
      + x0l5cff * dl5
      + x0l4cff * dl4
      + x1l3cff * y1 * dl13
      + x1l4cff * y1 * dl14
      + y1l3cff * y1 * y1 * dl13
      + y1l4cff * y1 * y1 * dl14;

    // The selected approximations for nf = 3,...,6
    double p3psapp1 = 0;
    double p3psapp2 = 0;
    if (nf == 3)
      {
        p3psapp1 = p3ps01
                   + 67731. * y1 * dl * ym
                   + 274100. * y1 * ym
                   - 104493. * y1 * (1.+2. * y)
                   + 34403. * y1 * y * y
                   + 353656. * y1 * dl
                   + 10620. * dl2
                   + 40006. * dl3
                   - 7412.1 * y1 * dl1
                   - 2365.1 * y1 * dl12
                   + 1533.0 * y1 * y1 * dl12;
        p3psapp2  = p3ps01
                    + 54593. * y1 * dl * ym
                    + 179748. * y1 * ym
                    - 195263. * y1
                    + 12789. * y1 * y * (1.+y)
                    + 4700.0 * y1 * dl
                    - 103604. * dl2
                    - 2758.3 * dl3
                    - 2801.2 * y1 * dl1
                    - 1986.9 * y1 * dl12
                    - 6005.9 * y1 * y1 * dl12;
      }
    else if (nf == 4)
      {
        p3psapp1 = p3ps01
                   + 90154. * y1 * dl * ym
                   + 359084. * y1 * ym
                   - 136319. * y1 * (1.+2. * y)
                   + 45379. * y1 * y * y
                   + 461167. * y1 * dl
                   + 13869. * dl2
                   + 52525. * dl3
                   - 7498.2 * y1 * dl1
                   - 2491.5 * y1 * dl12
                   + 1727.2 * y1 * y1 * dl12;
        p3psapp2  = p3ps01
                    + 72987. * y1 * dl * ym
                    + 235802. * y1 * ym
                    - 254921. * y1
                    + 17138. * y1 * y * (1.+y)
                    + 5212.9 * y1 * dl
                    - 135378. * dl2
                    - 3350.9 * dl3
                    - 1472.7 * y1 * dl1
                    - 1997.2 * y1 * dl12
                    - 8123.3 * y1 * y1 * dl12;
      }
    else if (nf == 5)
      {
        p3psapp1 = p3ps01
                   + 112481. * y1 * dl * ym
                   + 440555. * y1 * ym
                   - 166581. * y1 * (1.+2. * y)
                   + 56087. * y1 * y * y
                   + 562992. * y1 * dl
                   + 16882. * dl2
                   + 64577. * dl3
                   - 6570.1 * y1 * dl1
                   - 2365.7 * y1 * dl12
                   + 1761.7 * y1 * y1 * dl12;
        p3psapp2  = p3ps01
                    + 91468. * y1 * dl * ym
                    + 289658. * y1 * ym
                    - 311749. * y1
                    + 21521. * y1 * y * (1.+y)
                    + 4908.9 * y1 * dl
                    - 165795. * dl2
                    - 3814.9 * dl3
                    +  804.5 * y1 * dl1
                    - 1760.8 * y1 * dl12
                    - 10295. * y1 * y1 * dl12;
      }
    else if (nf == 6)
      {
        p3psapp1 = p3ps01
                   + 134701. * y1 * dl * ym
                   + 518318. * y1 * ym
                   - 195241. * y1 * (1.+2. * y)
                   + 66517. * y1 * y * y
                   + 658832. * y1 * dl
                   + 19605. * dl2
                   + 76125. * dl3
                   - 4734.5 * y1 * dl1
                   - 2035.2 * y1 * dl12
                   + 1633.1 * y1 * y1 * dl12;
        p3psapp2  = p3ps01
                    + 110032. * y1 * dl * ym
                    + 341158. * y1 * ym
                    - 365676. * y1
                    + 25934. * y1 * y * (1.+y)
                    + 3614.4 * y1 * dl
                    - 194868. * dl2
                    - 4172.2 * dl3
                    + 3924.3 * y1 * dl1
                    - 1324.9 * y1 * dl12
                    - 12520. * y1 * y1 * dl12;
      }
    else
      throw std::runtime_error(error("P3ps::Regular", "nf out of range."));

    // We return one of the two error-band representatives or their average
    if (_imod == 1)
      return p3psapp1;
    else if (_imod == 2)
      return p3psapp2;
    else
      return 0.5 * ( p3psapp1 + p3psapp2 );
  }

  // //_________________________________________________________________________________
  // P3qg::P3qg(int const& nf, double const& rho):
  //   Expression(),
  //   _nf(nf),
  //   _rho(rho)
  // {
  //   // Moments for the known exact small-x contribution (Vogt)
  //   const std::vector<double> N{-0.32821791 - _rho * 0.9999999999999999, -0.0079068 - _rho * 0.11111111111076043,
  //                               0.00668085 - _rho * 0.04000000000013782, 0.0112976 - _rho * 0.02040816};

  //   // Matrix
  //   const std::vector<std::vector<double>>  inv_A
  //   {
  //     {1.50859345e1, -1.73502577e2, 3.81181701e2, -2.25021001e2},
  //     {1.39350535e1, -2.33166358e2, 5.53133581e2, -3.37644014e2},
  //     {3.92339377, -8.52278788e1, 2.57493287e2, -1.79628628e2},
  //     {-1.15513419e-2, 2.74882781e-1, -9.19800075e-1, 7.20490480e-1}
  //   };

  //   // Matrix multiplication of inv_A and N
  //   _C.resize(N.size(), 0.);
  //   for (int i = 0; i < (int) N.size(); i++)
  //     for (int j = 0; j < (int) N.size(); j++)
  //       _C[i] += inv_A[i][j] * N[j];
  // }
  // double P3qg::Regular(double const& x) const
  // {
  //   const double dl  = log(x);
  //   const double dl2 = dl * dl;
  //   const double dlm = log(1 - x);
  //   return _nf * ( _C[0] * dl2 + _C[1] * dl + _C[2] * pow(x, 2) + _C[3] * pow(dlm, 4) - _rho * dl / x
  //                  + pow(CA, 3) * ( 82. / 81. + 2 * zeta3 ) * dl2 / x / 6 / pow(M_PI, 4) ) * pow(FourPi, 4);
  // }

  //_________________________________________________________________________________
  P3qg::P3qg(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double P3qg::Regular(double const& y) const
  {
    const int nf      = _nf;
    const int nf2     = _nf * _nf;
    const int nf3     = _nf * nf2;
    const double ym   = 1 / y;
    const double y1   = 1 - y;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double dl6  = dl * dl5;
    const double dl1  = log(1 - y);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;

    // Known large-x coefficients
    const double x1l5cff =   1.8518519e0 * nf - 4.1152263e-1 * nf2;
    const double x1l4cff =   3.5687794e1 * nf - 3.5116598e0 * nf2 - 8.2304527e-2 * nf3;
    const double y1l5cff =   2.8806584e0 * nf + 8.2304527e-1 * nf2;
    const double y1l4cff = - 4.0511391e1 * nf + 5.5418381e0 * nf2 + 1.6460905e-1 * nf3;

    // Known small-x coefficients
    const double bfkl1   =   3.9357613e3 * nf;
    const double x0l6cff = - 1.9588477e1 * nf + 2.7654321e0 * nf2;
    const double x0l5cff =   2.1573663e1 * nf + 1.7244444e1 * nf2;
    const double x0l4cff = - 2.8667643e3 * nf + 3.0122403e2 * nf2 + 4.1316872e0 * nf3 ;

    // The resulting part of the function
    const double p3qg01 =
      + bfkl1 * ym * dl2
      + x0l6cff * dl6
      + x0l5cff * dl5
      + x0l4cff * dl4
      + x1l4cff * dl14
      + x1l5cff * dl15
      + y1l4cff * y1 * dl14
      + y1l5cff * y1 * dl15;

    // The selected approximations for nf = 3,...,6
    double p3qgapp1 = 0;
    double p3qgapp2 = 0;
    if (nf == 3)
      {
        p3qgapp1 = p3qg01
                   + 187500. * ym * dl
                   + 826060. * ym * y1
                   - 150474.
                   + 226254. * y * (2.-y)
                   + 577733. * dl
                   - 180747. * dl2
                   + 95411. * dl3
                   +  119.8 * dl13
                   + 7156.3 * dl12
                   + 45790. * dl1
                   - 95682. * dl * dl1;
        p3qgapp2  = p3qg01
                    + 135000. * ym * dl
                    + 484742. * ym * y1
                    - 11627.
                    - 187478. * y * (2.-y)
                    + 413512. * dl
                    - 82500. * dl2
                    + 29987. * dl3
                    -  850.1 * dl13
                    - 11425. * dl12
                    - 75323. * dl1
                    + 282836. * dl * dl1;
      }
    else if (nf == 4)
      {
        p3qgapp1 = p3qg01
                   + 250000. * ym * dl
                   + 1089180. * ym * y1
                   - 241088.
                   + 342902. * y * (2.-y)
                   + 720081. * dl
                   - 247071. * dl2
                   + 126405. * dl3
                   +  272.4 * dl13
                   + 10911. * dl12
                   + 60563. * dl1
                   - 161448. * dl * dl1;
        p3qgapp2  = p3qg01
                    + 180000. * ym * dl
                    + 634090. * ym * y1
                    - 55958.
                    - 208744. * y * (2.-y)
                    + 501120. * dl
                    - 116073. * dl2
                    + 39173. * dl3
                    - 1020.8 * dl13
                    - 13864. * dl12
                    - 100922. * dl1
                    + 343243. * dl * dl1;
      }
    else if (nf == 5)
      {
        p3qgapp1 = p3qg01
                   + 312500. * ym * dl
                   + 1345700. * ym * y1
                   - 350466.
                   + 480028. * y * (2.-y)
                   + 837903. * dl
                   - 315928. * dl2
                   + 157086. * dl3
                   +  472.7 * dl13
                   + 15415. * dl12
                   + 75644. * dl1
                   - 244869. * dl * dl1;
        p3qgapp2  = p3qg01
                    + 225000. * ym * dl
                    + 776837. * ym * y1
                    - 119054.
                    - 209530. * y * (2.-y)
                    + 564202. * dl
                    - 152181. * dl2
                    + 48046. * dl3
                    - 1143.8 * dl13
                    - 15553. * dl12
                    - 126212. * dl1
                    + 385995. * dl * dl1;
      }
    else if (nf == 6)
      {
        p3qgapp1 = p3qg01
                   + 375000. * ym * dl
                   + 1595330. * ym * y1
                   - 477729.
                   + 637552. * y * (2.-y)
                   + 931556. * dl
                   - 387017. * dl2
                   + 187509. * dl3
                   +  715.5 * dl13
                   + 20710. * dl12
                   + 91373. * dl1
                   - 346374. * dl * dl1;
        p3qgapp2  = p3qg01
                    + 270000. * ym * dl
                    + 912695. * ym * y1
                    - 200034.
                    - 189918. * y * (2.-y)
                    + 603114. * dl
                    - 190521. * dl2
                    + 56661. * dl3
                    - 1224.3 * dl13
                    - 16453. * dl12
                    - 150856. * dl1
                    + 410661. * dl * dl1;
      }
    else
      throw std::runtime_error(error("P3qg::Regular", "nf out of range."));

    // We return one of the two error-band representatives or their average
    if (_imod == 1)
      return p3qgapp1;
    else if (_imod == 2)
      return p3qgapp2;
    else
      return 0.5 * ( p3qgapp1 + p3qgapp2 );
  }

  // //_________________________________________________________________________________
  // P3gq::P3gq(double const& rho):
  //   Expression(),
  //   _rho(rho)
  // {
  //   // Moments for the known exact small-x contribution (Vogt)
  //   const std::vector<double> N{-0.78602265 - _rho * 1.9999999999999991, 0.05776039 - _rho * 0.07407407407348075,
  //                               0.04868481 - _rho * 0.016000000001040217, 0.04071129 - _rho * 0.0058309037900882415};

  //   // Matrix
  //   const std::vector<std::vector<double>>  inv_A
  //   {
  //     {-15.48054909, -39.22467343, 13.63767859, 3.69946174},
  //     {442.68554891, 1038.70260338, -415.55730182,-113.79603728},
  //     {-1263.91174116,  -2938.1594667, 1302.32965951,362.67726047},
  //     {874.63704072, 2025.64293423,-958.33374045, -272.56485302}
  //   };

  //   // Matrix multiplication of inv_A and N
  //   _C.resize(N.size(), 0.);
  //   for (int i = 0; i < (int) N.size(); i++)
  //     for (int j = 0; j < (int) N.size(); j++)
  //       _C[i] += inv_A[j][i] * N[j];
  // }
  // double P3gq::Regular(double const& x) const
  // {
  //   const double dl  = log(x);
  //   const double dl2 = dl * dl;
  //   const double dl3 = dl * dl2;
  //   const double dlm = log(1 - x);
  //   return ( - _C[0] * dl / x + _C[1] * dl3 + _C[2] * x + _C[3] * dlm + _rho * dl2 / x
  //            - pow(CA, 3) * CF * zeta3 * dl3 / x / 3 / pow(M_PI, 4) ) * pow(FourPi, 4);
  // }

  //_________________________________________________________________________________
  P3gq::P3gq(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double P3gq::Regular(double const& y) const
  {
    const int nf      = _nf;
    const int nf2     = _nf * _nf;
    const double ym   = 1 / y;
    const double y1   = 1 - y;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double dl6  = dl * dl5;
    const double dl1  = log(1 - y);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;
    const double dl15 = dl1 * dl14;

// Known large-x coefficients
    const double x1l5cff = 1.3443073e+1 - 5.4869684e-1 * nf;
    const double x1l4cff = 3.7539831e+2 - 3.4494742e+1 * nf + 8.7791495e-1 * nf2;
    const double y1l5cff = 2.2222222e+1 - 5.4869684e-1 * nf;
    const double y1l4cff = 6.6242163e+2 - 4.7992684e+1 * nf + 8.7791495e-1 * nf2;

// x^-1 small-x coeffs, casimir scaled from p_gg (approx. for bfkl1)
    const double bfkl0 =   - 8.3086173e+3 / 2.25;
    const double bfkl1 = ( - 1.0691199e+5 - nf * 9.9638304e+2 ) / 2.25;

// Small-x double-logs with x^0
    const double x0l6cff =   5.2235940e+1 - 7.3744856e+0 * nf;
    const double x0l5cff = - 2.9221399e+2 + 1.8436214e+0 * nf;
    const double x0l4cff =   7.3106077e+3 - 3.7887135e+2 * nf - 3.2438957e+1 * nf2;

// The resulting part of the function
    const double p3gq01 =
      + bfkl0 * ym * dl3
      + bfkl1 * ym * dl2
      + x0l6cff * dl6
      + x0l5cff * dl5
      + x0l4cff * dl4
      + x1l4cff * dl14
      + x1l5cff * dl15
      + y1l4cff * y1 * dl14
      + y1l5cff * y1 * dl15;

// The selected approximations for nf = 3,...,6
    double p3gqapp1 = 0;
    double p3gqapp2 = 0;
    if (nf == 3)
      {
        p3gqapp1 = p3gq01
                   + 3.5 * bfkl1 * ym * dl
                   - 27891. * ym * y1
                   - 309124.
                   + 1056866. * y * (2.-y)
                   - 124735. * dl
                   - 16246. * dl2
                   + 131175. * dl3
                   + 4970.1 * dl13
                   + 60041. * dl12
                   + 343181. * dl1
                   - 958330. * dl * dl1;
        p3gqapp2 = p3gq01
                   + 7. * bfkl1 * ym * dl
                   - 1139334. * ym * y1
                   + 143008.
                   - 290390. * y * (2.-y)
                   - 659492. * dl
                   + 303685. * dl2
                   - 81867. * dl3
                   + 1811.8 * dl13
                   -  465.9 * dl12
                   - 51206. * dl1
                   + 274249. * dl * dl1;
      }
    else if (nf == 4)
      {
        p3gqapp1 = p3gq01
                   + 3.5 * bfkl1 * ym * dl
                   - 8302.8 * ym * y1
                   - 347706.
                   + 1105306. * y * (2.-y)
                   - 127650. * dl
                   - 29728. * dl2
                   + 137537. * dl3
                   + 4658.1 * dl13
                   + 59205. * dl12
                   + 345513. * dl1
                   - 995120. * dl * dl1;
        p3gqapp2 = p3gq01
                   + 7. * bfkl1 * ym * dl
                   - 1129822. * ym * y1
                   + 108527.
                   - 254166. * y * (2.-y)
                   - 667254. * dl
                   + 293099. * dl2
                   - 77437. * dl3
                   + 1471.3 * dl13
                   - 1850.3 * dl12
                   - 52451. * dl1
                   + 248634. * dl * dl1;
      }
    else if (nf == 5)
      {
        p3gqapp1 = p3gq01
                   + 3.5 * bfkl1 * ym * dl
                   + 14035. * ym * y1
                   - 384003.
                   + 1152711. * y * (2.-y)
                   - 126346. * dl
                   - 42967. * dl2
                   + 144270. * dl3
                   + 4385.5 * dl13
                   + 58688. * dl12
                   + 348988. * dl1
                   - 1031165. * dl * dl1;
        p3gqapp2  = p3gq01
                    + 7. * bfkl1 * ym * dl
                    - 1117561. * ym * y1
                    + 76329.
                    - 218973. * y * (2.-y)
                    - 670799. * dl
                    + 282763. * dl2
                    - 72633. * dl3
                    + 1170.0 * dl13
                    - 2915.5 * dl12
                    - 52548. * dl1
                    + 223771. * dl * dl1;
      }
    else if (nf == 6)
      {
        p3gqapp1 = p3gq01
                   + 3.5 * bfkl1 * ym * dl
                   + 39203. * ym * y1
                   - 417914.
                   + 1199042. * y * (2.-y)
                   - 120750. * dl
                   - 55941. * dl2
                   + 151383. * dl3
                   + 4149.2 * dl13
                   + 58466. * dl12
                   + 353589. * dl1
                   - 1066510. * dl * dl1;
        p3gqapp2  = p3gq01
                    + 7. * bfkl1 * ym * dl
                    - 1102470. * ym * y1
                    + 46517.
                    - 184858. * y * (2.-y)
                    - 670056. * dl
                    + 272689. * dl2
                    - 67453. * dl3
                    + 905.0 * dl13
                    - 3686.2 * dl12
                    - 51523. * dl1
                    + 199594. * dl * dl1;
      }
    else
      throw std::runtime_error(error("P3gq::Regular", "nf out of range."));

// We return one of the two error-band representatives or their average
    if (_imod == 1)
      return p3gqapp1;
    else if (_imod == 2)
      return p3gqapp2;
    else
      return 0.5 * ( p3gqapp1 + p3gqapp2 );
  }

  // //_________________________________________________________________________________
  // P3gg::P3gg(double const& rho):
  //   Expression(),
  //   _rho(rho)
  // {
  //   // Moments for the known exact small-x contribution (Vogt)
  //   const std::vector<double> N{6.94542399 - _rho * 0.9999999999999999, 0.01255781 - _rho * 0.11111111111076043,
  //                               -0.24160401 - _rho * 0.04000000000013782, -0.27416992 - _rho * 0.02040816};

  //   // Matrix
  //   const std::vector<std::vector<double>>  inv_A
  //   {
  //     {1.57592805e1, -1.89525930e2, 4.34798303e2, -2.67019531e2},
  //     {1.51046749e1, -2.60999383e2, 6.46267166e2, -4.10596673e2},
  //     {6.05043429, -1.35844233e2, 4.26863395e2, -3.12298274e2},
  //     {-3.86434270e-1, 9.19582568, -3.07706474e1, 2.41030189e1}
  //   };

  //   // Matrix multiplication of inv_A and N
  //   _C.resize(N.size(), 0.);
  //   for (int i = 0; i < (int) N.size(); i++)
  //     for (int j = 0; j < (int) N.size(); j++)
  //       _C[i] += inv_A[i][j] * N[j];
  // }
  // double P3gg::Regular(double const& x) const
  // {
  //   const double nf  = 4;
  //   const double dl  = log(x);
  //   const double dl2 = dl * dl;
  //   const double dl3 = dl * dl2;
  //   const double dlm = log(1 - x);
  //   return ( _C[0] * dl2 + _C[1] * dl + _C[2] * pow(x, 2) + _C[3] * pow(dlm, 2) - _rho * dl / x
  //            - pow(CA, 4) / (3 * pow(Pi2, 2)) * zeta3 * dl3 / x
  //            + ( pow(CA, 4) * ( - 1205. / 162. + 67. / 36. * zeta2 + 1. / 4. * pow(zeta2, 2) - 11. / 2. * zeta3 )
  //                + nf * pow(CA, 3) * ( - 233. / 162. + 13. / 36. * zeta2 - 1. / 3. * zeta3 )
  //                + nf * pow(CA, 2) * CF * ( 617. / 243. - 13. / 18. * zeta2 + 2. / 3. * zeta3) ) * dl2 / x / pow(Pi2, 2) / 2 ) * pow(FourPi, 4);
  // }

  // //_________________________________________________________________________________
  // P3gg::P3gg(int const& nf, int const& imod):
  //   Expression(),
  //   _nf(std::min(nf, 5)),
  //   _imod(imod)
  // {
  //   _A4gluon = 40880.330e0 - 11714.246e0 * _nf + 440.04876e0 * pow(_nf, 2) + 7.3627750e0 * pow(_nf, 3);
  // }
  // double P3gg::Regular(double const& x) const
  // {
  //   const int nf2     = _nf * _nf;
  //   const double xm   = 1 / x;
  //   const double x1   = 1 - x;
  //   const double dl   = log(x);
  //   const double dl2  = dl * dl;
  //   const double dl3  = dl * dl2;
  //   const double dlm  = log(1 - x);
  //   const double dlm2 = dlm * dlm;
  //   const double dlm3 = dlm * dlm2;

  //   // The known large-x coefficients [except delta(1-x)]
  //   const double Ccoeff  = 8.5814120e4 - 1.3880515e4 * _nf + 1.3511111e2 * nf2;
  //   const double Dcoeff  = 5.4482808e4 - 4.3411337e3 * _nf - 2.1333333e1 * nf2;

  //   // The known coefficients of 1/x*ln^a x terms, a = 3,2
  //   const double bfkl0 = - 8.308617314e3;
  //   const double bfkl1 = - 1.069119905e5 - 9.963830436e2 * _nf;

  //   // The resulting part of the function
  //   const double P3gg01 =
  //     + bfkl0 * dl3 * xm
  //     + bfkl1 * dl2 * xm
  //     + Ccoeff * dlm
  //     + Dcoeff - _A4gluon;

  //   // The selected approximations for nf = 3, 4, 5
  //   double P3ggApp1 = P3gg01;
  //   double P3ggApp2 = P3gg01;
  //   if (_nf <= 3)
  //     {
  //       P3ggApp1 +=
  //         + 3.4 * bfkl1 * dl * xm
  //         - 345063. * x1 * xm
  //         + 86650. * ( 1 + x * x ) * x1
  //         + 158160. * dl
  //         - 15741. * x1 * dlm2
  //         - 9417. * x1 * dlm3;
  //       P3ggApp2 +=
  //         + 5.4 * bfkl1 * dl * xm
  //         - 1265632. * x1 * xm
  //         - 656644. * ( 1 + x * x ) * x1
  //         - 1352233. * dl
  //         + 203298. * x1 * dlm2
  //         + 39112. * x1 * dlm3;
  //     }
  //   else if (_nf == 4)
  //     {
  //       P3ggApp1 +=
  //         + 3.4 * bfkl1 * dl * xm
  //         - 342625. * x1 * xm
  //         + 100372. * ( 1 + x * x ) * x1
  //         + 189167. * dl
  //         - 29762. * x1 * dlm2
  //         - 12102. * x1 * dlm3;
  //       P3ggApp2 +=
  //         + 5.4 * bfkl1 * dl * xm
  //         - 1271540. * x1 * xm
  //         - 649661. * ( 1 + x * x ) * x1
  //         - 1334919. * dl
  //         + 191263. * x1 * dlm2
  //         + 36867. * x1 * dlm3;
  //     }
  //   else if (_nf >= 5)
  //     {
  //       P3ggApp1 +=
  //         + 3.4 * bfkl1 * dl * xm
  //         - 337540. * x1 * xm
  //         + 119366. * ( 1 + x * x ) * x1
  //         + 223769. * dl
  //         - 45129. * x1 * dlm2
  //         - 15046. * x1 * dlm3;
  //       P3ggApp2 +=
  //         + 5.4 * bfkl1 * dl * xm
  //         - 1274800. * x1 * xm
  //         - 637406. * ( 1 + x * x ) * x1
  //         - 1314010. * dl
  //         + 177882. * x1 * dlm2
  //         + 34362. * x1 * dlm3;
  //     }

  //   // We return (for now) one of the two error-band boundaries or the
  //   // present best estimate, their average
  //   if (_imod == 1)
  //     return P3ggApp1;
  //   else if (_imod == 2)
  //     return P3ggApp2;
  //   else
  //     return 0.5 * ( P3ggApp1 + P3ggApp2 );
  // }
  // double P3gg::Singular(double const& x) const
  // {
  //   return _A4gluon / ( 1 - x );
  // }
  // double P3gg::Local(double const& x) const
  // {
  //   const double B4gluon = 68587.64 - 18143.983e0 * _nf + 423.81135e0 * pow(_nf, 2) + 9.0672154e-1 * pow(_nf, 3);
  //   return log(1 - x) * _A4gluon + B4gluon;
  // }

  //_________________________________________________________________________________
  P3gg::P3gg(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
    _A4gluon = 40880.330e0 - 11714.246e0 * _nf + 440.04876e0 * pow(_nf, 2) + 7.3627750e0 * pow(_nf, 3);
  }
  double P3gg::Regular(double const& y) const
  {
    const int nf      = _nf;
    const int nf2     = _nf * _nf;
    const int nf3     = _nf * nf2;
    const double ym   = 1 / y;
    const double y1   = 1 - y;
    const double dl   = log(y);
    const double dl2  = dl * dl;
    const double dl3  = dl * dl2;
    const double dl4  = dl * dl3;
    const double dl5  = dl * dl4;
    const double dl6  = dl * dl5;
    const double dl1  = log(1 - y);
    const double dl12 = dl1 * dl1;
    const double dl13 = dl1 * dl12;
    const double dl14 = dl1 * dl13;

    // Known large-y coefficients [except delta(1-x)]
    const double ccoeff  = 8.5814120e4 - 1.3880515e4 * nf + 1.3511111e2 * nf2;
    const double dcoeff  = 5.4482808e4 - 4.3411337e3 * nf - 2.1333333e1 * nf2;

    const double x1l4cff = 5.6460905e1 * nf - 3.6213992e0 * nf2;
    const double x1l3cff = 2.4755054e2 * nf - 4.0559671e1 * nf2 + 1.5802469e0 * nf3;

    // Known small-x coefficients
    const double bfkl0  = - 8.3086173e3;
    const double bfkl1  = - 1.0691199e5 - 9.9638304e2 * nf;

    const double x0l6cff =  1.44e2 - 2.7786008e1 * nf + 7.9012346e-1 * nf2;
    const double x0l5cff = -1.44e2 - 1.6208066e2 * nf + 1.4380247e1 * nf2;
    const double x0l4cff =  2.6165784e4 - 3.3447551e3 * nf + 9.1522635e1 * nf2 - 1.9753086e-1 * nf3;

    // The resulting part of the function
    const double p3gg01 =
      + bfkl0 * dl3 * ym
      + bfkl1 * dl2 * ym
      + x0l6cff * dl6
      + x0l5cff * dl5
      + x0l4cff * dl4
      + ccoeff * dl1
      + dcoeff - _A4gluon
      + x1l4cff * y1 * dl14
      + x1l3cff * y1 * dl13;

    // The selected approximations for nf = 3,...,6
    double p3ggapp1 = 0;
    double p3ggapp2 = 0;
    if (nf == 3)
      {
        p3ggapp1 = p3gg01
                   - 421311. * y1 * dl * ym
                   - 325557. * y1 * ym
                   + 1679790. * y1
                   - 1456863. * y1 * y
                   + 3246307. * y1 * dl
                   + 2026324. * dl * dl
                   + 549188. * dl3
                   +   8337. * y1 * dl1
                   +  26718. * y1 * dl1 * dl1
                   -  27049. * y1 * y1 * dl13;
        p3ggapp2 = p3gg01
                   - 700113. * y1 * dl * ym
                   - 2300581. * y1 * ym
                   + 896407. * y1 * (1.+2. * y)
                   - 162733. * y1 * y * y
                   - 2661862. * y1 * dl
                   + 196759. * dl * dl
                   - 260607. * dl3
                   +  84068. * y1 * dl1
                   + 346318. * y1 * dl1 * dl1
                   + 315725. * dl * dl1 * dl1;
      }
    else if (nf == 4)
      {
        p3ggapp1 = p3gg01
                   - 437084. * y1 * dl * ym
                   - 361570. * y1 * ym
                   + 1696070. * y1
                   - 1457385. * y1 * y
                   + 3195104. * y1 * dl
                   + 2009021. * dl * dl
                   + 544380. * dl3
                   +  9938. * y1 * dl1
                   +  24376. * y1 * dl1 * dl1
                   -  22143. * y1 * y1 * dl13;
        p3ggapp2 = p3gg01
                   - 706649. * y1 * dl * ym
                   - 2274637. * y1 * ym
                   + 836544. * y1 * (1.+2. * y)
                   - 199929. * y1 * y * y
                   - 2683760. * y1 * dl
                   + 168802. * dl * dl
                   - 250799. * dl3
                   +  36967. * y1 * dl1
                   +  24530. * y1 * dl1 * dl1
                   -  71470. * y1 * y1 * dl1 * dl1;
      }
    else if (nf == 5)
      {
        p3ggapp1 = p3gg01
                   - 439426. * y1 * dl * ym
                   - 293679. * y1 * ym
                   + 1916281. * y1
                   - 1615883. * y1 * y
                   + 3648786. * y1 * dl
                   + 2166231. * dl * dl
                   + 594588. * dl3
                   +  50406. * y1 * dl1
                   +  24692. * y1 * dl1 * dl1
                   + 174067. * y1 * y1 * dl1;
        p3ggapp2 = p3gg01
                   - 705978. * y1 * dl * ym
                   - 2192234. * y1 * ym
                   + 1730508. * y1 * y
                   + 353143. * y1 * (2.-y * y)
                   - 2602682. * y1 * dl
                   + 178960. * dl * dl
                   - 218133. * dl3
                   +   2285. * y1 * dl1
                   +  19295. * y1 * dl1 * dl1
                   -  13719. * y1 * y1 * dl1 * dl1;
      }
    else if (nf == 6)
      {
        p3ggapp1 = p3gg01
                   - 476018. * y1 * dl * ym
                   - 469289. * y1 * ym
                   + 2049351. * y1
                   - 1589000. * y1 * y
                   + 3185549. * y1 * dl
                   + 1994521. * dl * dl
                   + 527723. * dl3
                   - 340674. * y1 * dl1
                   +  22460. * y1 * dl1 * dl1
                   - 394556. * dl * dl1;
        p3ggapp2 = p3gg01
                   - 709863. * y1 * dl * ym
                   - 2134347. * y1 * ym
                   + 1605315. * y1 * y
                   + 360743. * y1 * (2.-y * y)
                   - 2426250. * y1 * dl
                   + 230631. * dl * dl
                   - 185804. * dl3
                   - 7992.9 * y1 * dl1
                   + 15918. * y1 * dl1 * dl1
                   - 32771. * y1 * y1 * dl1;
      }
    else
      throw std::runtime_error(error("P3gg::Regular", "nf out of range."));

    // We return one of the two error-band representatives or their average
    if (_imod == 1)
      return p3ggapp1;
    else if (_imod == 2)
      return p3ggapp2;
    else
      return 0.5 * ( p3ggapp1 + p3ggapp2 );
  }
  double P3gg::Singular(double const& y) const
  {
    return _A4gluon / ( 1 - y );
  }
  double P3gg::Local(double const& y) const
  {
    double B4gluon = 68587.64 - 18143.983e0 * _nf + 423.81135e0 * pow(_nf, 2) + 9.0672154e-1 * pow(_nf, 3);
    if (_imod == 1)
      B4gluon += - 0.2;
    if (_imod == 2)
      B4gluon += 0.2;
    return log(1 - y) * _A4gluon + B4gluon;
  }
}
