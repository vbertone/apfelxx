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
  aP3nsp::aP3nsp(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double aP3nsp::Regular(double const& y) const
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
  double aP3nsp::Singular(double const& y) const
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
  double aP3nsp::Local(double const& y) const
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
  aP3nsm::aP3nsm(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double aP3nsm::Regular(double const& y) const
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
  double aP3nsm::Singular(double const& y) const
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
  double aP3nsm::Local(double const& y) const
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
  aP3nss::aP3nss(int const& nf, int const& imod):
    Expression(),
    _nf(nf),
    _imod(imod)
  {
  }
  double aP3nss::Regular(double const& y) const
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

  //_________________________________________________________________________________
  P3nsp::P3nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CF4    = CF * CF3;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double nf2    = _nf * _nf;
    const double nf3    = _nf * nf2;
    const double d4FF   = 5. / 36.;
    const double d4FA   = 5. / 2.;
    const double zeta32 =  zeta3 * zeta3;

    // Coefficient of the singular part
    _csing = nf3*CF*(-0.3950617283950617 + (64*zeta3)/27.) + nf2*CA*CF*(11.395061728395062 - (608*zeta2)/81. + (2240*zeta3)/27. - (112*zeta4)/3.) + CF2*nf2*(29.530864197530864 - (640*zeta3)/9. + 32*zeta4) + d4FF*nf*(256*zeta2 - (256*zeta3)/3. - (1280*zeta5)/3.) + CF3*nf*(63.55555555555556 + (592*zeta3)/3. - 320*zeta5) + CF2*CA*nf*(-420.5679012345679 + (440*zeta2)/3. + (3712*zeta3)/9. - 128*zeta2*zeta3 - 176*zeta4 + 160*zeta5) + CA2*CF*nf*(-297.98765432098764 + (20320*zeta2)/81. - (23104*zeta3)/27. + (448*zeta2*zeta3)/3. - (176*zeta4)/3. + (2096*zeta5)/9.) + d4FA*(-384*zeta32 - 128*zeta2 + (128*zeta3)/3. + (3520*zeta5)/3. - 992*zeta6) + CA3*CF*(1040.469135802469 - 16*zeta32 - (88400*zeta2)/81. + (20944*zeta3)/27. - (352*zeta2*zeta3)/3. + 1804*zeta4 - (3608*zeta5)/9. - (2504*zeta6)/3.);

    // Coefficient of the local part
    _cloc = nf3*CF*(-1.617283950617284 + (32*zeta2)/81. + (304*zeta3)/81. - (32*zeta4)/27.) + nf2*CA*CF*(-3.574074074074074 + (3170*zeta2)/81. - (320*zeta3)/9. + (80*zeta2*zeta3)/3. - (80*zeta4)/9. - (88*zeta5)/9.) + CF2*nf2*(-6.962962962962963 + (1244*zeta2)/27. + (56*zeta3)/27. - (160*zeta2*zeta3)/9. - (2104*zeta4)/27. + (368*zeta5)/9.) + CF3*nf*(32 + 224*zeta32 + 162*zeta2 - 308*zeta3 - (256*zeta2*zeta3)/3. - 204*zeta4 + 912*zeta5 - (6434*zeta6)/9.) + CA2*CF*nf*(185.4351851851852 + (416*zeta32)/3. - (41092*zeta2)/81. + (9554*zeta3)/27. - (580*zeta2*zeta3)/3. + (1234*zeta4)/9. + (1130*zeta5)/9. - (3913*zeta6)/27.) + d4FF*nf*(-192 + 256*zeta32 + (1888*zeta2)/3. - (992*zeta3)/3. + 64*zeta2*zeta3 - (352*zeta4)/3. - 1120*zeta5 + (1792*zeta6)/9.) + CF2*CA*nf*(-143.53703703703704 - (1232*zeta32)/3. - (3892*zeta2)/27. - (15400*zeta3)/27. + (2672*zeta2*zeta3)/9. + (27854*zeta4)/27. - (7432*zeta5)/9. + 351*zeta6) + CF3*CA*(-521.25 + 3220*zeta32 + 1167*zeta2 - 3260*zeta3 - (1988*zeta2*zeta3)/3. + 2167*zeta4 + 128*zeta3*zeta4 - 976*zeta5 + 2064*zeta2*zeta5 + (79297*zeta6)/18. - 10920*zeta7) + CA3*CF*(-576.841049382716 + (1672*zeta32)/3. + (13864*zeta2)/9. - (152284*zeta3)/81. + (584*zeta2*zeta3)/3. - (11206*zeta4)/27. + (32*zeta3*zeta4)/3. + (11522*zeta5)/9. + (1472*zeta2*zeta5)/3. + (73333*zeta6)/108. - (8960*zeta7)/3.) + d4FA*(96 - 704*zeta32 - (944*zeta2)/3. - (1232*zeta3)/3. - 896*zeta2*zeta3 + (32*zeta4)/3. - 64*zeta3*zeta4 - 400*zeta5 + 320*zeta2*zeta5 - (1562*zeta6)/9. + 2800*zeta7) + CF4*(203.04166666666666 - 1152*zeta32 - 450*zeta2 + 2004*zeta3 - 120*zeta2*zeta3 - 342*zeta4 + 64*zeta3*zeta4 - 2520*zeta5 - 384*zeta2*zeta5 - 2111*zeta6 + 5880*zeta7) + CA2*CF2*(823.3055555555555 - (7102*zeta32)/3. - (46771*zeta2)/27. + (129662*zeta3)/27. + (2096*zeta2*zeta3)/9. - (60850*zeta4)/27. - 32*zeta3*zeta4 + (5354*zeta5)/9. - 2104*zeta2*zeta5 - (5497*zeta6)/2. + 8610*zeta7);
  }
  double P3nsp::Regular(double const& x) const
  {
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CF4    = CF * CF3;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double nf2    = _nf * _nf;
    const double nf3    = _nf * nf2;
    const double d4FF   = 5. / 36.;
    const double d4FA   = 5. / 2.;
    const double x2     = x * x;
    const double x3     = x * x2;
    const double x4     = x * x3;
    const double x5     = x * x4;
    const double x6     = x * x5;
    const double x7     = x * x6;
    const double x8     = x * x7;
    const double omx    = 1 - x;
    const double omx2   = omx * omx;
    const double omx3   = omx * omx2;
    const double dl     = log(x);
    const double dl2    = dl * dl;
    const double dl3    = dl * dl2;
    const double dl4    = dl * dl3;
    const double dl5    = dl * dl4;
    const double dl6    = dl * dl5;
    const double dlm    = log(omx);
    const double dlm2   = dlm * dlm;
    const double dlm3   = dlm * dlm2;
    const double zeta32 = zeta3 * zeta3;
    return CF4*(121.73518255221838 - 137761.66937038003*x4 + 888660.0564771394*x5 - 255257.20831827258*x6 + 56504.70015666814*x7 - 7278.481839298749*x8 + x3*(-196836.46298324046 + 1.0996580796252204e7*dl2 - 6.798217475425384e6*dl3 + 1.340842882472595e6*dl4 - 229248.56457283554*dl5 + 12605.764064403846*dl6 - 2.0837808397693438e6*dl) + x2*(-180757.92110024643 - 1.1427983746250909e7*dl2 - 5.0079164792098645e6*dl3 - 1.058688186339151e6*dl4 - 99625.29984005733*dl5 - 6009.869754559037*dl6 + 1.5050014971501704e6*dl) + omx2*(16586.159495222895*dlm2 + 1038.401448631465*dlm3 + 110402.60213600783*dlm) + omx3*(-13972.072122209684*dlm2 + 14939.171721197328*dlm3 + 273767.4277486546*dlm) + (-397.1540209993292*dlm2 + 0.24353717515070297*dlm3 + 547.3217280230208*dlm)*omx + (-187793.43566668688 - 272975.4083456475*dl2 - 32222.89268614284*dl3 - 1899.367309507996*dl4 - 65.60390271679701*dl5 - 1.4220980561652938*dl6 - 992616.5325918229*dl)*x) + CA2*CF2*(-728.668509396418 - 167997.72967007212*x4 + 986146.2332668182*x5 - 268462.9724140581*x6 + 57771.64073654107*x7 - 7240.601318488816*x8 + x3*(-220193.98305849134 + 8.824375823232185e6*dl2 - 5.7398593715630155e6*dl3 + 1.0189101632009936e6*dl4 - 194250.00358779266*dl5 + 7112.390846749811*dl6 - 2.0242373129066303e6*dl) + x2*(-202657.24740802383 - 9.58711278625441e6*dl2 - 4.210101667003363e6*dl3 - 907810.8681686072*dl4 - 86078.89413077169*dl5 - 5323.493376729169*dl6 + 1.1890735285366746e6*dl) + omx2*(12769.367663618052*dlm2 + 768.3119553170851*dlm3 + 88450.7605898829*dlm) + omx3*(-10494.56996911999*dlm2 + 12172.812302881772*dlm3 + 222673.13555037434*dlm) + (-661.8195813012476*dlm2 + 0.18526411962171732*dlm3 - 1828.180282374631*dlm)*omx + (-200222.11550904167 - 250087.51217235738*dl2 - 28924.23478513794*dl3 - 1662.5749133856066*dl4 - 59.564825403137895*dl5 - 1.8185973174100227*dl6 - 911280.7829857245*dl)*x) + CF2*CA*_nf*(370.6001594432628 + 799434.5974829213*x4 - 67907.67129934975*x5 + 9880.545611107753*x6 - 1269.9942552438415*x7 + 74.68927728188892*x8 + x2*(-99174.85084372171 - 2.9898023369751773e6*dl2 - 631426.1187814261*dl3 - 58592.83992113247*dl4 - 3625.9742502930985*dl5 - 7.2653375741749e6*dl) + x3*(-141788.3896555702 - 3.808543091966822e6*dl2 + 777618.9527348245*dl3 - 126281.62912098368*dl4 + 8558.917188265248*dl5 + 5.634684668728132e6*dl) + omx2*(570.9521875047908*dlm2 + 35.307625663134054*dlm3 + 3095.7651560163004*dlm) + omx3*(-513.3493365021654*dlm2 + 455.29211536712467*dlm3 + 8333.830827089663*dlm) + (56.33453802583154*dlm2 + 0.009608500861648456*dlm3 + 695.0138554664134*dlm)*omx + (-495541.4430229262 - 19889.30196603619*dl2 - 1278.464112494632*dl3 - 50.8156065570392*dl4 - 2.7841754100636877*dl5 - 155251.58696595568*dl)*x) + d4FA*(-506.66340818310766 - 4.329556468716501e6*x4 + 508981.3459170827*x5 - 117868.25604538301*x6 + 24606.770421637448*x7 - 3074.2575871542062*x8 + x2*(1.9398652900259148e6 + 636585.9649376993*dl2 + 78148.53843811863*dl3 + 5009.325260994501*dl4 + 2.1318103492321437e6*dl) + x3*(1.8875659942023659e6 + 1.9415560919244129e6*dl2 + 197398.6790655809*dl3 + 78444.8677876684*dl4 + 3.576973766379004e6*dl) + omx2*(5171.07848405067*dlm2 + 331.7446262231697*dlm3 + 35250.80135330049*dlm) + omx3*(-4036.3487792030483*dlm2 + 4677.461593764502*dlm3 + 87539.14104842396*dlm) + (75.10476849634433*dlm2 + 0.13044891131346173*dlm3 + 607.1944877801232*dlm)*omx + (84132.64443824427 - 1877.13988474203*dl2 - 38.17797656984607*dl3 + 14.510990721390725*dl4 - 3351.563779565061*dl)*x) + nf2*CA*CF*(107.92551889087873 - 297.70764883478756*x4 + 42.01399152473712*x5 - 11.596855466174036*x6 + 2.7534590854193453*x7 - 0.4025822296503464*x8 + x2*(72.83488400209473 + 21.530781944043902*dl2 + 1.1809166879530217*dl3 + 0.08973776498862204*dl4 + 129.28699197826893*dl) + x3*(67.78094670252176 + 69.38350967404082*dl2 + 22.694030912388943*dl3 + 3.163423879086563*dl4 + 357.24999308537286*dl) + omx2*(1.8217092380978652*dlm2 + 0.11272803608622221*dlm3 + 14.233747014427246*dlm) + omx3*(-1.4171756023213926*dlm2 + 1.6241938240599578*dlm3 + 28.236966075945485*dlm) + (0.0009608074399317006*dlm2 + 0.000027234310987381596*dlm3 - 9.469801733439*dlm)*omx + (-8.089603268815807 - 1.6733275892250057*dl2 + 3.8034174084181855*dl3 + 0.44446537065492814*dl4 - 41.04322718063796*dl)*x) + nf3*CF*(2.849320066876777 + 1.5778955408839035*x4 - 0.2428128949159311*x5 + 0.07946833215589176*x6 - 0.021458020712737814*x7 + 0.003981190661424573*x8 + x2*(7.240534816215909 - 0.716558816419646*dl2 - 0.1811756547208233*dl3 + 2.275202396134914*dl) + x3*(-7.954453420016102 - 1.5838160725737989*dl2 + 0.03281266693551782*dl3 + 3.124356372847891*dl) + omx3*(0.057419921457276014*dlm2 - 0.04004790832050069*dlm3 - 0.7450079808395396*dlm) + omx2*(-0.05022057865224671*dlm2 - 0.0031700366899497004*dlm3 - 0.3392957609314987*dlm) + (-0.00002953650315556271*dlm2 - 8.434884097324452e-7*dlm3 - 0.00035650148080391205*dlm)*omx + (-1.7127538494132848 + 0.09888998684399462*dl2 - 0.0987619258909032*dl3 + 2.766957117784918*dl)*x) + CF2*nf2*(-31.327270381911053 + 557.0249813832176*x4 - 79.30115418703222*x5 + 22.01533518002141*x6 - 5.23966601700404*x7 + 0.7625913053893202*x8 + x3*(-174.99499169459878 - 148.50059831055054*dl2 - 34.93164233841145*dl3 - 7.285176288384152*dl4 - 648.1631918432682*dl) + x2*(-184.34582034713034 - 20.308939110350888*dl2 + 9.894262353740865*dl3 + 1.418247640083144*dl4 - 227.2960956698377*dl) + omx3*(3.6325667862153983*dlm2 - 2.9718169465160837*dlm3 - 52.105746356010826*dlm) + omx2*(-3.4575561146559917*dlm2 - 0.21419487155635558*dlm3 + 12.207229364996438*dlm) + (-0.0018263228124021519*dlm2 - 0.0000517557761867042*dlm3 - 27.676526507241057*dlm)*omx + (-259.95126953131216 - 5.201696577727609*dl2 - 5.8266275240497345*dl3 - 0.4444391803746076*dl4 + 13.872433818871743*dl)*x) + d4FF*_nf*(-123.8949105144963 + 7756.561787567016*x4 - 1465.8144670310546*x5 + 397.4749710700201*x6 - 89.6714067694069*x7 + 10.829318643634043*x8 + x3*(15768.109243027382 + 1359.6345015828995*dl2 - 897.6356307396716*dl3 - 19453.035062524057*dl) + x2*(-22730.285780160026 - 776.1597099178529*dl2 - 57.75829747940975*dl3 - 9132.81086646499*dl) + omx3*(49.30497145936295*dlm2 + 15.061672795339389*dlm3 + 98.97222878925778*dlm) + omx2*(2.4434834042739535*dlm2 - 0.5242569579866798*dlm3 + 143.53928915894295*dlm) + (103.38634822188118*dlm2 - 0.0021996168803625166*dlm3 - 550.5998005914264*dlm)*omx + (86.21392136313702 - 0.38558578640141494*dl2 - 0.010760648217204617*dl3 + 1586.1598305321193*dl)*x) + CA2*CF*_nf*(-887.3724270672252 - 204987.45583928737*x4 + 16901.44195032383*x5 - 2378.6683590547*x6 + 289.6530398882654*x7 - 14.28663874844101*x8 + x3*(36442.01740471579 + 974023.537118974*dl2 - 198046.03785629736*dl3 + 32265.218772401502*dl4 - 2154.146609880292*dl5 - 1.439892304434672e6*dl) + x2*(25516.511338875236 + 765273.7259452463*dl2 + 161941.7424526509*dl3 + 15035.704235225456*dl4 + 933.1658219000408*dl5 + 1.859548642988758e6*dl) + omx3*(137.4696382374695*dlm2 - 126.46582822197561*dlm3 - 2322.638153793919*dlm) + omx2*(-152.30666318928016*dlm2 - 9.53264165810479*dlm3 - 1024.5106646046727*dlm) + (-6.442431739430879*dlm2 - 0.0024636803429999394*dlm3 - 134.05338833196546*dlm)*omx + (128897.32027199857 + 5155.70717501893*dl2 + 328.46130521796755*dl3 + 11.240309110575737*dl4 + 0.7021449010525752*dl5 + 40048.15200000848*dl)*x) + CF3*_nf*(113.41221457835701 - 767091.3027906149*x4 + 67353.01196871749*x5 - 10058.271981251883*x6 + 1339.072385984154*x7 - 84.76953943607374*x8 + x3*(138577.39123400202 + 3.6823290441043833e6*dl2 - 756023.1375171966*dl3 + 122202.07382561707*dl4 - 8439.911701988003*dl5 - 5.459801343849226e6*dl) + x2*(97446.95746981385 + 2.8862623960465887e6*dl2 + 607955.7080382529*dl3 + 56341.6979506201*dl4 + 3472.729891374825*dl5 + 7.011348675681651e6*dl) + omx3*(717.606547128552*dlm2 - 518.7784148927017*dlm3 - 9800.22111124406*dlm) + omx2*(-677.025573480704*dlm2 - 45.20002399972181*dlm3 - 4645.676081926284*dlm) + (-155.16913512094587*dlm2 - 0.01472552973385878*dlm3 - 306.6155867170534*dlm)*omx + (468632.1777954987 + 19090.085881841096*dl2 + 1242.2641338231685*dl3 + 49.341578928705054*dl4 + 2.305824070339687*dl5 + 148619.10906279602*dl)*x) + CA3*CF*(1369.216887096375 + 34118.92218901031*x4 - 211373.19275504685*x5 + 57734.501139348504*x6 - 12398.489671128014*x7 + 1540.3997104452503*x8 + x2*(42567.74198550853 + 2.3095183359807394e6*dl2 + 1.0105275555364704e6*dl3 + 217078.7587967714*dl4 + 20521.31958156479*dl5 + 1263.7164471829185*dl6 - 291911.46194368193*dl) + x3*(46351.00517614079 - 2.159531725563034e6*dl2 + 1.3712314052620581e6*dl3 - 253429.3997975875*dl4 + 46221.71835164751*dl5 - 1971.7482860179384*dl6 + 465496.8109809639*dl) + omx3*(1724.5119193371845*dlm2 - 2210.3009088962867*dlm3 - 40390.87947077632*dlm) + omx2*(-2205.9873120443294*dlm2 - 134.09881258930523*dlm3 - 15304.527003547117*dlm) + (161.1867148029785*dlm2 - 0.03461957826159862*dlm3 + 424.1546800696331*dlm)*omx + (43012.48613470983 + 57719.618832880726*dl2 + 6704.8563777494655*dl3 + 386.80597907635683*dl4 + 14.982073423927492*dl5 + 0.48375285002626467*dl6 + 207289.66976113827*dl)*x) + CF3*CA*(-692.1934770585552 + 265696.88491901703*x4 - 1.5968128439986897e6*x5 + 442983.58123761095*x6 - 96367.17450121965*x7 + 12220.624733525234*x8 + x2*(327962.37316182157 + 1.6962648731186297e7*dl2 + 7.445803409061612e6*dl3 + 1.5931792133384976e6*dl4 + 150680.08194620957*dl5 + 9228.434334895213*dl6 - 2.152293192045573e6*dl) + x3*(356494.3161569876 - 1.5871225912189448e7*dl2 + 1.0142445565402087e7*dl3 - 1.8727925363546691e6*dl4 + 343015.05029526824*dl5 - 14889.673374939224*dl6 + 3.401481149806218e6*dl) + omx3*(20629.888349439283*dlm2 - 22758.686382907068*dlm3 - 416181.2637399787*dlm) + omx2*(-24521.261272747648*dlm2 - 1469.8523879596457*dlm3 - 165095.6880263144*dlm) + (1338.0618050064645*dlm2 - 0.34857985087752774*dlm3 + 868.3501374635733*dlm)*omx + (328675.4849485309 + 429472.2963143897*dl2 + 50033.50514405304*dl3 + 2919.268383126514*dl4 + 100.44156542235578*dl5 + 2.5619362962574574*dl6 + 1.566757506902313e6*dl)*x) + nf3*CF*(-0.3950617283950617 - (88*dl2)/81. - (8*dl3)/81. - (64*dl)/27. - (32*zeta3)/9.) + CF2*nf2*(16.691358024691358 + (152*dl3)/27. + (4*dl4)/9. + (1216*dlm)/81. + (3040*zeta2)/81. + dl2*(26.469135802469136 + (16*zeta2)/9.) + dl*(95.50617283950618 + (64*zeta2)/9.) + (3296*zeta3)/27. - (320*zeta4)/9.) + nf2*CA*CF*(31.69753086419753 + (44*dl3)/27. + dl2*(23.85185185185185 - (16*zeta2)/3.) - (632*zeta2)/81. - (3616*zeta3)/27. + dl*(93.34567901234568 - (304*zeta2)/9. + (32*zeta3)/3.) + (488*zeta4)/9.) + CF2*CA*_nf*(-451.41975308641975 - (44*dl4)/9. + dl3*(-77.48148148148148 - 48*zeta2) - (117152*zeta2)/81. + dlm*(-512.395061728395 + (640*zeta2)/3. - (896*zeta3)/3.) + dl2*(-407.01234567901236 - (4928*zeta2)/9. - 64*zeta3) - (45712*zeta3)/27. + (128*zeta2*zeta3)/3. + dl*(-1117.7530864197531 - (9584*zeta2)/9. - (1328*zeta3)/3. - (1448*zeta4)/3.) + (5272*zeta4)/9. - (872*zeta5)/3.) + CF3*_nf*(-72 - (20*dl4)/9. - (4*dl5)/9. + (5272*zeta2)/3. + dl3*(-5.333333333333333 + 64*zeta2) + (1600*zeta3)/3. + (256*zeta2*zeta3)/3. + dl2*(-56.666666666666664 + (1840*zeta2)/3. + (256*zeta3)/3.) + dlm*(-293.3333333333333 + 256*zeta3) - 488*zeta4 + dl*(166.66666666666666 + (3592*zeta2)/3. + (2080*zeta3)/3. + 304*zeta4) + (1616*zeta5)/3.) + CF3*CA*(-1119.6666666666667 + (22*dl5)/9. - 1928*zeta32 + (129704*zeta2)/3. + dl4*(18.88888888888889 + 168*zeta2) - 4828*zeta3 - (53200*zeta2*zeta3)/3. + dl3*(3.3333333333333335 + 240*zeta2 + 48*zeta3) + dl2*(98 + (3032*zeta2)/3. + (416*zeta3)/3. - 7972*zeta4) - 20594*zeta4 + dl*(-920.3333333333334 + (47312*zeta2)/3. - (12448*zeta3)/3. - 13248*zeta2*zeta3 - 16976*zeta4 - 960*zeta5) - 3556*zeta5 - 39158*zeta6) + CA3*CF*(266.7037037037037 + (128*zeta32)/3. + (1363420*zeta2)/81. + 38*dl4*zeta2 + dl3*(16.432098765432098 + (320*zeta2)/3.) - (90592*zeta3)/27. - 4792*zeta2*zeta3 + dl*(1808.4197530864199 + (37580*zeta2)/9. - 264*zeta3 - 3280*zeta2*zeta3 - 3368*zeta4) + dl2*(308.679012345679 + (866*zeta2)/3. - 1040*zeta4) - (89947*zeta4)/9. - (2090*zeta5)/9. - (45661*zeta6)/6.) + d4FF*_nf*(128 - 512*zeta32 + (3904*zeta2)/3. + (4480*zeta3)/3. - 256*zeta2*zeta3 - 672*zeta4 + (1280*zeta5)/3. - 1984*zeta6) + CA2*CF*_nf*(-177.56172839506172 - (64*zeta32)/3. - (116*zeta2)/81. + dl3*(-8.962962962962964 + 16*zeta2) + dl2*(-154.44444444444446 + (620*zeta2)/3.) + (43520*zeta3)/27. - (544*zeta2*zeta3)/3. + (1028*zeta4)/9. + dl*(-796.0617283950618 + (6104*zeta2)/9. - (32*zeta3)/3. + 56*zeta4) - (2600*zeta5)/9. - (248*zeta6)/3.) + CF4*(10.666666666666666 - (4*dl5)/3. + dl6/9. + 1984*zeta32 + dl4*(5.333333333333333 - (248*zeta2)/3.) - (37172*zeta2)/3. + dl3*(-35.333333333333336 - 64*zeta2 - (320*zeta3)/3.) - (2240*zeta3)/3. + 5792*zeta2*zeta3 + 10588*zeta4 + dl2*(-112 - 1528*zeta2 - 256*zeta3 + 5320*zeta4) + (4040*zeta5)/3. + dl*(-130 - 6200*zeta2 - 944*zeta3 + 6144*zeta2*zeta3 + 8624*zeta4 + 1920*zeta5) + (42376*zeta6)/3.) + d4FA*(-64 + 1568*zeta32 - (70720*zeta2)/3. - 256*dl3*zeta2 - 48*dl4*zeta2 - (9376*zeta3)/3. + 9664*zeta2*zeta3 + 11200*zeta4 + dl2*(-2112*zeta2 + 1920*zeta4) + dl*(-9984*zeta2 + 6144*zeta2*zeta3 + 5760*zeta4) - (560*zeta5)/3. + 15388*zeta6) + CA2*CF2*(3298.283950617284 + 540*zeta32 + dl3*(241.85185185185185 - 360*zeta2) + dl4*(13.444444444444445 - 140*zeta2) - (4007276*zeta2)/81. + (184004*zeta3)/27. + (49984*zeta2*zeta3)/3. + (192076*zeta4)/9. + dlm*(2193.382716049383 - (4288*zeta2)/3. + (704*zeta3)/3. + 864*zeta4) + dl2*(1416.5432098765432 - (12404*zeta2)/9. + 32*zeta3 + 4812*zeta4) + dl*(3138.58024691358 - (149132*zeta2)/9. + (8984*zeta3)/3. + 11232*zeta2*zeta3 + (40388*zeta4)/3. - 240*zeta5) + (8396*zeta5)/3. + (101669*zeta6)/3.);
  }
  double P3nsp::Singular(double const& x) const
  {
    return _csing / ( 1 - x );
  }
  double P3nsp::Local(double const& x) const
  {
    return _cloc + _csing * log(1 - x);
  }

  //_________________________________________________________________________________
  P3nsm::P3nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CF4    = CF * CF3;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double nf2    = _nf * _nf;
    const double nf3    = _nf * nf2;
    const double d4FF   = 5. / 36.;
    const double d4FA   = 5. / 2.;
    const double zeta32 = zeta3 * zeta3;

    // Coefficient of the singular part
    _csing = nf3*CF*(-0.3950617283950617 + (64*zeta3)/27.) + nf2*CA*CF*(11.395061728395062 - (608*zeta2)/81. + (2240*zeta3)/27. - (112*zeta4)/3.) + CF2*nf2*(29.530864197530864 - (640*zeta3)/9. + 32*zeta4) + d4FF*nf*(256*zeta2 - (256*zeta3)/3. - (1280*zeta5)/3.) + CF3*nf*(63.55555555555556 + (592*zeta3)/3. - 320*zeta5) + CF2*CA*nf*(-420.5679012345679 + (440*zeta2)/3. + (3712*zeta3)/9. - 128*zeta2*zeta3 - 176*zeta4 + 160*zeta5) + CA2*CF*nf*(-297.98765432098764 + (20320*zeta2)/81. - (23104*zeta3)/27. + (448*zeta2*zeta3)/3. - (176*zeta4)/3. + (2096*zeta5)/9.) + d4FA*(-384*zeta32 - 128*zeta2 + (128*zeta3)/3. + (3520*zeta5)/3. - 992*zeta6) + CA3*CF*(1040.469135802469 - 16*zeta32 - (88400*zeta2)/81. + (20944*zeta3)/27. - (352*zeta2*zeta3)/3. + 1804*zeta4 - (3608*zeta5)/9. - (2504*zeta6)/3.);

    // Coefficient of the local part
    _cloc = nf3*CF*(-1.617283950617284 + (32*zeta2)/81. + (304*zeta3)/81. - (32*zeta4)/27.) + nf2*CA*CF*(-3.574074074074074 + (3170*zeta2)/81. - (320*zeta3)/9. + (80*zeta2*zeta3)/3. - (80*zeta4)/9. - (88*zeta5)/9.) + CF2*nf2*(-6.962962962962963 + (1244*zeta2)/27. + (56*zeta3)/27. - (160*zeta2*zeta3)/9. - (2104*zeta4)/27. + (368*zeta5)/9.) + CF3*nf*(32 + 224*zeta32 + 162*zeta2 - 308*zeta3 - (256*zeta2*zeta3)/3. - 204*zeta4 + 912*zeta5 - (6434*zeta6)/9.) + CA2*CF*nf*(185.4351851851852 + (416*zeta32)/3. - (41092*zeta2)/81. + (9554*zeta3)/27. - (580*zeta2*zeta3)/3. + (1234*zeta4)/9. + (1130*zeta5)/9. - (3913*zeta6)/27.) + d4FF*nf*(-192 + 256*zeta32 + (1888*zeta2)/3. - (992*zeta3)/3. + 64*zeta2*zeta3 - (352*zeta4)/3. - 1120*zeta5 + (1792*zeta6)/9.) + CF2*CA*nf*(-143.53703703703704 - (1232*zeta32)/3. - (3892*zeta2)/27. - (15400*zeta3)/27. + (2672*zeta2*zeta3)/9. + (27854*zeta4)/27. - (7432*zeta5)/9. + 351*zeta6) + CF3*CA*(-521.25 + 3220*zeta32 + 1167*zeta2 - 3260*zeta3 - (1988*zeta2*zeta3)/3. + 2167*zeta4 + 128*zeta3*zeta4 - 976*zeta5 + 2064*zeta2*zeta5 + (79297*zeta6)/18. - 10920*zeta7) + CA3*CF*(-576.841049382716 + (1672*zeta32)/3. + (13864*zeta2)/9. - (152284*zeta3)/81. + (584*zeta2*zeta3)/3. - (11206*zeta4)/27. + (32*zeta3*zeta4)/3. + (11522*zeta5)/9. + (1472*zeta2*zeta5)/3. + (73333*zeta6)/108. - (8960*zeta7)/3.) + d4FA*(96 - 704*zeta32 - (944*zeta2)/3. - (1232*zeta3)/3. - 896*zeta2*zeta3 + (32*zeta4)/3. - 64*zeta3*zeta4 - 400*zeta5 + 320*zeta2*zeta5 - (1562*zeta6)/9. + 2800*zeta7) + CF4*(203.04166666666666 - 1152*zeta32 - 450*zeta2 + 2004*zeta3 - 120*zeta2*zeta3 - 342*zeta4 + 64*zeta3*zeta4 - 2520*zeta5 - 384*zeta2*zeta5 - 2111*zeta6 + 5880*zeta7) + CA2*CF2*(823.3055555555555 - (7102*zeta32)/3. - (46771*zeta2)/27. + (129662*zeta3)/27. + (2096*zeta2*zeta3)/9. - (60850*zeta4)/27. - 32*zeta3*zeta4 + (5354*zeta5)/9. - 2104*zeta2*zeta5 - (5497*zeta6)/2. + 8610*zeta7);
  }
  double P3nsm::Regular(double const& x) const
  {
    const double CF2    = CF * CF;
    const double CF3    = CF * CF2;
    const double CF4    = CF * CF3;
    const double CA2    = CA * CA;
    const double CA3    = CA * CA2;
    const double nf2    = _nf * _nf;
    const double nf3    = _nf * nf2;
    const double d4FF   = 5. / 36.;
    const double d4FA   = 5. / 2.;
    const double x2     = x * x;
    const double x3     = x * x2;
    const double x4     = x * x3;
    const double x5     = x * x4;
    const double x6     = x * x5;
    const double x7     = x * x6;
    const double x8     = x * x7;
    const double omx    = 1 - x;
    const double omx2   = omx * omx;
    const double omx3   = omx * omx2;
    const double dl     = log(x);
    const double dl2    = dl * dl;
    const double dl3    = dl * dl2;
    const double dl4    = dl * dl3;
    const double dl5    = dl * dl4;
    const double dl6    = dl * dl5;
    const double dlm    = log(omx);
    const double dlm2   = dlm * dlm;
    const double dlm3   = dlm * dlm2;
    const double zeta32 = zeta3 * zeta3;
    return CF3*_nf*(113.4120011174599 + 787667.5006855791*x4 - 68944.41842103818*x5 + 10444.757676331326*x6 - 1429.3517136370274*x7 + 100.90309206931782*x8 + x2*(-96519.37592063358 - 2.9340805919559924e6*dl2 - 618512.3410903295*dl3 - 57330.94696879176*dl4 - 3538.5528631916163*dl5 - 7.135628003351676e6*dl) + x3*(-138350.49163423278 - 3.740043929391008e6*dl2 + 766122.8569936334*dl3 - 124147.99008644585*dl4 + 8520.171490437902*dl5 + 5.529008686694059e6*dl) + omx2*(266.5359487842269*dlm2 + 13.867860498918636*dlm3 + 1776.9187575857466*dlm) + omx3*(-190.03579309525333*dlm2 + 252.05562965820678*dlm3 + 4475.263419548655*dlm) + (-154.64370308875067*dlm2 + 0.0002093661917441519*dlm3 - 354.97007323696033*dlm)*omx + (-492194.9725145989 - 19141.699732073655*dl2 - 1230.4312227300568*dl3 - 41.506554179371335*dl4 - 1.0765269645845665*dl5 - 152426.3024361366*dl)*x) + CA3*CF*(1369.217519745186 - 3.167881711882702e6*x4 + 348119.3683051468*x5 - 77368.55928350042*x6 + 15688.824444186066*x7 - 1929.8497157624467*x8 + x2*(1.439824098101666e6 + 526299.661927352*dl2 + 68481.40052896808*dl3 + 4371.5079671261465*dl4 + 1.6836547425229268e6*dl) + x3*(1.4027632858755412e6 + 1.5195805473395456e6*dl2 + 135221.43730526615*dl3 + 63070.51393958626*dl4 + 2.512513110795898e6*dl) + omx2*(4177.310114904137*dlm2 + 252.80594247052406*dlm3 + 28856.27811679158*dlm) + omx3*(-3887.409565287275*dlm2 + 3888.607346724233*dlm3 + 70972.98833072656*dlm) + (164.1565305097854*dlm2 + 0.048774245433928595*dlm3 + 609.795490950153*dlm)*omx + (36225.29391643096 - 3169.0079796313844*dl2 + 137.9995460570473*dl3 + 61.46292492602557*dl4 - 25284.961198493354*dl)*x) + CF2*CA*_nf*(370.6003037488185 - 462508.2416989307*x4 + 24447.28902853565*x5 + 330.9263809726774*x6 - 1048.4823118340987*x7 + 281.1092423149506*x8 + x3*(259838.6172544616 + 365330.62775616406*dl2 + 5483.820939105746*dl3 + 15785.155341676613*dl4 + 174513.055220927*dl) + x2*(263279.26213149785 + 133657.25041702559*dl2 + 18962.27321648177*dl3 + 1360.0550938462907*dl4 + 342100.3265687574*dl) + omx3*(3184.3110552344992*dlm2 - 2570.082689755375*dlm3 - 47518.097787198836*dlm) + omx2*(-3066.5582174804717*dlm2 - 192.0036436947431*dlm3 - 21629.536475379886*dlm) + (54.31702504129046*dlm2 - 0.04775744473275153*dlm3 + 791.2961359608357*dlm)*omx + (-87758.01293199288 - 988.0495485093652*dl2 - 82.06315662638144*dl3 - 5.278253811335637*dl4 - 6107.9967223525755*dl)*x) + d4FF*_nf*(-123.89491025696954 - 11531.986154326696*x4 + 1994.4421549201586*x5 - 600.9716183870293*x6 + 148.25775779875013*x7 - 23.207216688958674*x8 + x2*(29136.96598793167 + 3213.0393361286797*dl2 + 337.8877003544379*dl3 + 14382.228098150168*dl) + x3*(-20350.66858316156 - 1213.7485250357531*dl2 + 1178.8950976332733*dl3 + 26474.205228356233*dl) + omx2*(125.75298965964336*dlm2 + 7.049420659444683*dlm3 + 989.8725879460665*dlm) + omx3*(-59.821788972978645*dlm2 + 126.03013723657844*dlm3 + 2124.460713664157*dlm) + (103.44825997834246*dlm2 - 0.00045180719369195953*dlm3 - 549.8443183952782*dlm)*omx + (1552.8175713892379 - 0.26010203052752384*dl2 - 0.0066212493037045754*dl3 - 112.97222573861521*dl)*x) + nf2*CA*CF*(107.92551891691646 + 345.32988032993137*x4 - 52.28200320120174*x5 + 15.686683844082285*x6 - 3.95008763455563*x7 + 0.64263966466164*x8 + x3*(1594.2194527773343 + 170.57224874226287*dl2 - 57.323720836592976*dl3 - 1381.5099070572926*dl) + x2*(-1788.292719111005 - 164.6052941973352*dl2 - 10.537025904433872*dl3 - 932.9200650027723*dl) + omx3*(6.25729838749635*dlm2 - 4.396691644027627*dlm3 - 82.76591642264498*dlm) + omx2*(-5.445464178167316*dlm2 - 0.3424019447490028*dlm3 - 35.09204367873701*dlm) + (-0.003132685552463405*dlm2 - 0.00008932225315034865*dlm3 - 5.9637934245955035*dlm)*omx + (-367.442354199891 - 4.440986860249375*dl2 + 1.6305741303485624*dl3 - 68.50311472577103*dl)*x) + nf3*CF*(2.849320066876777 + 1.5778955408839035*x4 - 0.2428128949159311*x5 + 0.07946833215589176*x6 - 0.021458020712737814*x7 + 0.003981190661424573*x8 + x2*(7.240534816215909 - 0.716558816419646*dl2 - 0.1811756547208233*dl3 + 2.275202396134914*dl) + x3*(-7.954453420016102 - 1.5838160725737989*dl2 + 0.03281266693551782*dl3 + 3.124356372847891*dl) + omx3*(0.057419921457276014*dlm2 - 0.04004790832050069*dlm3 - 0.7450079808395396*dlm) + omx2*(-0.05022057865224671*dlm2 - 0.0031700366899497004*dlm3 - 0.3392957609314987*dlm) + (-0.00002953650315556271*dlm2 - 8.434884097324452e-7*dlm3 - 0.00035650148080391205*dlm)*omx + (-1.7127538494132848 + 0.09888998684399462*dl2 - 0.0987619258909032*dl3 + 2.766957117784918*dl)*x) + CF2*nf2*(-31.327270397119623 - 481.08894631542955*x4 + 79.52756210961634*x5 - 22.85735731471352*x6 + 5.5065485801989995*x7 - 0.8136351975238172*x8 + x2*(122.93305400257658 + 2.1883445849374894*dl2 - 1.4926248019419401*dl3 - 0.11039609496027437*dl4 + 139.68588899894712*dl) + x3*(113.42212556293273 + 36.403606452595085*dl2 + 49.55995032597273*dl3 - 0.05660354544904601*dl4 + 657.038315653278*dl) + omx2*(3.723226793364631*dlm2 + 0.23088701511423912*dlm3 + 61.23937734284941*dlm) + omx3*(-3.1721502487369735*dlm2 + 3.2631949654841494*dlm3 + 61.68644504185927*dlm) + (0.001983598289117164*dlm2 + 0.00005625570897168679*dlm3 - 34.741331774717764*dlm)*omx + (306.6908527139879 + 0.46258752550226667*dl2 - 1.4759460881820659*dl3 + 0.44455634898561014*dl4 + 70.21398376725462*dl)*x) + CA2*CF*_nf*(-887.3724602781684 - 27169.426501192134*x4 + 7991.108241202866*x5 - 4062.7896168319885*x6 + 1286.492888078763*x7 - 254.2687314315045*x8 + x2*(-178712.12913621156 - 2733.1887734703046*dl2 - 234.56547158478818*dl3 - 16081.744208728664*dl) + x3*(135543.99809670684 + 6455.691916855717*dl2 - 1297.2910736715903*dl3 - 12991.615329976765*dl) + omx2*(2977.9306742459303*dlm2 + 187.37815357365395*dlm3 + 20164.8076716133*dlm) + omx3*(-3304.801645755466*dlm2 + 2406.766935131687*dlm3 + 44670.275995742675*dlm) + (-4.637323182222187*dlm2 + 0.04901212005553582*dlm3 - 158.92219830571142*dlm)*omx + (68288.64669366983 + 204.94184313672926*dl2 + 17.240188957721895*dl3 + 890.2537335813006*dl)*x) + d4FA*(-506.6635780077066 + 4.622114690252967e6*x4 - 538439.8919028032*x5 + 126716.46836564715*x6 - 26923.278918400705*x7 + 3528.0126001461804*x8 + x3*(-1.969515379221794e6 - 2.0164382703023518e6*dl2 - 220693.35578531132*dl3 - 82944.98012873487*dl4 - 3.9475154768777634e6*dl) + x2*(-2.0262887772176547e6 - 691887.5489957804*dl2 - 89272.25992109484*dl3 - 5824.5249867833345*dl4 - 2.3330426765472745e6*dl) + omx3*(11023.015966757292*dlm2 - 10567.28110793159*dlm3 - 192045.35092650057*dlm) + omx2*(-11483.66526330298*dlm2 - 687.6110596362346*dlm3 - 79312.6513020516*dlm) + (66.89212025151248*dlm2 - 0.1011662205949918*dlm3 + 506.8795273387555*dlm)*omx + (-197133.66207966537 + 956.9933151215305*dl2 - 351.3282831987833*dl3 - 77.61354982878737*dl4 + 45754.07924075656*dl)*x) + CA2*CF2*(-728.6712270831409 + 1.4241417283367464e7*x4 - 1.5685079256841217e6*x5 + 348765.8835480323*x6 - 70736.07381296491*x7 + 8686.83251256517*x8 + x3*(-6.306106568756023e6 - 6.830482178074538e6*dl2 - 607233.0686799908*dl3 - 282905.1456679518*dl4 - 1.1283478629923292e7*dl) + x2*(-6.472861368971234e6 - 2.3605279385382356e6*dl2 - 305759.81975707575*dl3 - 19510.803193158154*dl4 - 7.5593930570589565e6*dl) + omx3*(15895.27169451532*dlm2 - 16597.769752564756*dlm3 - 302593.4066956257*dlm) + omx2*(-17513.37831457797*dlm2 - 1069.8147361125834*dlm3 - 120818.27283989845*dlm) + (-676.0326162292066*dlm2 - 0.2141140159883218*dlm3 - 2451.314411759972*dlm)*omx + (-174191.0723327886 + 13399.355607032627*dl2 - 446.0023807975091*dl3 - 212.20681109865998*dl4 + 93251.76729806514*dl)*x) + CF3*CA*(-692.1893416387634 - 9.962241394257875e6*x4 + 2.2920324243658893e6*x5 - 566465.1556926225*x6 + 120267.634449529*x7 - 15450.955922185132*x8 + x3*(2.360959760440094e6 + 5.53076064866208e7*dl2 - 1.364408905960994e7*dl3 + 1.9056761631784623e6*dl4 - 222961.44237587837*dl5 - 8.566014616781591e7*dl) + x2*(1.7739047933254014e6 + 4.160260431364982e7*dl2 + 7.8969469695762e6*dl3 + 707550.1555457156*dl4 + 36443.42183609345*dl5 + 1.0125747827535492e8*dl) + omx2*(43062.903355551585*dlm2 + 2668.8590936989794*dlm3 + 299580.1569490727*dlm) + omx3*(-41208.04899037176*dlm2 + 39108.61570781123*dlm3 + 716874.1094195545*dlm) + (1371.5386664651849*dlm2 + 0.5959962169983869*dlm3 + 1578.063952545771*dlm)*omx + (3.9926481299965256e6 + 98464.40379486867*dl2 + 7691.925041041126*dl3 + 472.7735343480742*dl4 + 5.4002689994334725*dl5 + 852618.969499146*dl)*x) + CF4*(121.73047345037378 + 154747.88277554477*x4 - 956192.7576718889*x5 + 270300.05198926176*x6 - 59756.33041132081*x7 + 7769.780726916616*x8 + x2*(190032.33160055583 + 9.483117480970714e6*dl2 + 4.1926739586004033e6*dl3 + 892937.4646463047*dl4 + 84555.44861057342*dl5 + 5142.049126439339*dl6 - 1.1877773418158703e6*dl) + x3*(206612.29691312896 - 8.860216992281096e6*dl2 + 5.745848068750033e6*dl3 - 1.047518912736338e6*dl4 + 195072.02668046625*dl5 - 8350.779124625957*dl6 + 1.954415823507961e6*dl) + omx3*(18606.037687600412*dlm2 - 18382.58397326322*dlm3 - 336497.6625129542*dlm) + omx2*(-19715.57127121846*dlm2 - 1183.2296251357725*dlm3 - 139343.02810175167*dlm) + (-415.0657521561032*dlm2 - 0.26170700033084665*dlm3 + 328.5699195257726*dlm)*omx + (188166.56945165712 + 252131.1006306902*dl2 + 29405.253696824544*dl3 + 1801.4508313381757*dl4 + 62.68043791283*dl5 + 1.0016401872997682*dl6 + 939491.6489762529*dl)*x) + nf3*CF*(-0.3950617283950617 - (88*dl2)/81. - (8*dl3)/81. - (64*dl)/27. - (32*zeta3)/9.) + CF2*nf2*(-92.74074074074075 - (664*dl3)/81. - (4*dl4)/9. + (1216*dlm)/81. + dl2*(-32.592592592592595 - (16*zeta2)/3.) + (1952*zeta2)/81. + dl*(-23.80246913580247 - (1024*zeta2)/27. - (64*zeta3)/9.) + (608*zeta3)/27. - 32*zeta4) + nf2*CA*CF*(86.41358024691358 + (692*dl3)/81. + (4*dl4)/9. + dl2*(53.382716049382715 - (16*zeta2)/9.) - (88*zeta2)/81. - (2272*zeta3)/27. + dl*(153 - (304*zeta2)/27. + (128*zeta3)/9.) + (472*zeta4)/9.) + CF2*CA*_nf*(-3578.925925925926 - (68*dl4)/3. - (32*dl5)/15. + (51800*zeta2)/81. + dl3*(-52.79012345679013 + (368*zeta2)/9.) + dlm*(-512.395061728395 + (640*zeta2)/3. - (896*zeta3)/3.) + (77888*zeta3)/27. + (1568*zeta2*zeta3)/3. + dl2*(-507.7037037037037 + 592*zeta2 + (704*zeta3)/3.) + (3416*zeta4)/3. + dl*(-3107.8765432098767 + (42416*zeta2)/27. + (18032*zeta3)/9. + 224*zeta4) - (2344*zeta5)/3.) + CF3*_nf*(2904 + (196*dl4)/9. + (76*dl5)/45. + dl3*(119.11111111111111 - (128*zeta2)/3.) - 904*zeta2 + dl2*(895.3333333333334 - (1744*zeta2)/3. - (320*zeta3)/3.) - (5504*zeta3)/3. - (832*zeta2*zeta3)/3. + dlm*(-293.3333333333333 + 256*zeta3) + dl*(2918.6666666666665 - (2936*zeta2)/3. - 1248*zeta3 - 416*zeta4) - 104*zeta4 + 432*zeta5) + CA2*CF2*(3125.8888888888887 + (104*dl5)/15. - (8*dl6)/9. - 4284*zeta32 + (1011020*zeta2)/81. + dl4*(98.11111111111111 + 132*zeta2) + (89012*zeta3)/27. - (38192*zeta2*zeta3)/3. + dl3*(-49.01234567901235 + (2776*zeta2)/9. + (1280*zeta3)/3.) + dl2*(274.962962962963 - (2468*zeta2)/3. + (4960*zeta3)/3. - 2172*zeta4) + (11296*zeta4)/3. + dlm*(2193.382716049383 - (4288*zeta2)/3. + (704*zeta3)/3. + 864*zeta4) + (18172*zeta5)/3. + dl*(6693.3456790123455 + (86372*zeta2)/27. - (9976*zeta3)/9. - 7536*zeta2*zeta3 - 2712*zeta4 + 4096*zeta5) - (34277*zeta6)/3.) + CF4*(-12816 + (92*dl5)/15. - (19*dl6)/45. - 3008*zeta32 + (55652*zeta2)/3. + dl4*(18.666666666666668 + (376*zeta2)/3.) + 12392*zeta3 - 8320*zeta2*zeta3 + dl3*(183.33333333333334 + 256*zeta2 + (1088*zeta3)/3.) + dl2*(1024 + 2040*zeta2 + 768*zeta3 - 3592*zeta4) - 1816*zeta4 - 2616*zeta5 + dl*(552.6666666666666 + (27848*zeta2)/3. + 5152*zeta3 - 7168*zeta2*zeta3 - (21688*zeta4)/3. + 1760*zeta5) - 8888*zeta6) + d4FF*_nf*(-3776 - 512*zeta32 + (3584*zeta2)/3. + 2112*zeta3 - 256*zeta2*zeta3 + (3136*zeta4)/3. + dl*(-1322.6666666666667 + (640*zeta2)/3. + 384*zeta3 + (1024*zeta4)/3.) + 1280*zeta5 - 1984*zeta6) + d4FA*(949.3333333333334 + 928*zeta32 - (6752*zeta2)/3. + 16*dl4*zeta2 + 4480*zeta3 - 2688*zeta2*zeta3 + dl2*(-64*zeta2 - 640*zeta4) + (38608*zeta4)/3. + dl*(192 - 640*zeta2 + 2112*zeta3 - 384*zeta2*zeta3 + 1360*zeta4 - 4000*zeta5) - 4496*zeta5 - (2860*zeta6)/3.) + CA2*CF*_nf*(723.5246913580247 + (26*dl4)/9. + (8*dl5)/15. - (64*zeta32)/3. + dl3*(-52.41975308641975 - (16*zeta2)/9.) - (30520*zeta2)/81. + dl2*(-342.0987654320988 - (580*zeta2)/9. - (304*zeta3)/3.) - (2644*zeta3)/27. - (992*zeta2*zeta3)/3. + dl*(-461.44444444444446 - (2704*zeta2)/27. - (6808*zeta3)/9. - (1120*zeta4)/9.) - 294*zeta4 - (104*zeta5)/3. - (248*zeta6)/3.) + CA3*CF*(1734.0123456790122 - (8*dl5)/15. + (4*dl6)/15. + (3256*zeta32)/3. + dl4*(-9.333333333333334 - 30*zeta2) - (121276*zeta2)/81. + dl3*(235.64197530864197 - (992*zeta2)/9. - (256*zeta3)/3.) - (34264*zeta3)/27. + 3376*zeta2*zeta3 - (42961*zeta4)/9. + dl2*(1556.8024691358025 - (1270*zeta2)/9. - (1136*zeta3)/3. + (1388*zeta4)/3.) + dl*(2329.703703703704 - (33652*zeta2)/27. + (9200*zeta3)/9. + 1456*zeta2*zeta3 - (58*zeta4)/3. - (3580*zeta5)/3.) - (8078*zeta5)/3. + (67057*zeta6)/18.) + CF3*CA*(-399.6666666666667 - (586*dl5)/45. + (44*dl6)/45. + 6152*zeta32 + dl4*(-119.77777777777777 - (656*zeta2)/3.) - 26380*zeta2 + dl3*(-401.1111111111111 - 432*zeta2 - (2096*zeta3)/3.) - (41920*zeta3)/3. + (52480*zeta2*zeta3)/3. - 400*zeta4 + dl2*(-3179.3333333333335 - (1520*zeta2)/3. - (6304*zeta3)/3. + 4868*zeta4) + dl*(-10488.333333333334 - (34136*zeta2)/3. - 4488*zeta3 + 13088*zeta2*zeta3 + (31916*zeta4)/3. - 4112*zeta5) + 2468*zeta5 + (61082*zeta6)/3.);
  }
  double P3nsm::Singular(double const& x) const
  {
    return _csing / ( 1 - x );
  }
  double P3nsm::Local(double const& x) const
  {
    return _cloc + _csing * log(1 - x);
  }

  //_________________________________________________________________________________
  P3nss::P3nss(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P3nss::Regular(double const& x) const
  {
    const double nf2    = _nf * _nf;
    const double d33c   = 5. / 18.;
    const double zeta32 = zeta3 * zeta3;
    const double x2     = x * x;
    const double x3     = x * x2;
    const double x4     = x * x3;
    const double x5     = x * x4;
    const double x6     = x * x5;
    const double x7     = x * x6;
    const double x8     = x * x7;
    const double omx    = 1 - x;
    const double omx2   = omx * omx;
    const double omx3   = omx * omx2;
    const double dl     = log(x);
    const double dl2    = dl * dl;
    const double dl3    = dl * dl2;
    const double dl4    = dl * dl3;
    const double dl5    = dl * dl4;
    const double dl6    = dl * dl5;
    const double dlm    = log(omx);
    const double dlm2   = dlm * dlm;
    const double dlm3   = dlm * dlm2;
    return nf2*d33c*(2.2372847169819517e-7 + 1597.1970598101527*x4 - 196.00864203198256*x5 + 43.831055972943524*x6 - 9.277810393710654*x7 + 1.3461215825361428*x8 + x3*(10751.684725885321 + 853.0162380265076*dl2 - 305.83143843874046*dl3 - 7519.694114104166*dl) + x2*(-12782.959507592881 - 1299.394376974993*dl2 - 1.234542687938335*dl3 - 5641.709178519317*dl) + omx2*(18.073466366678986*dlm2 + 1.0031692988720067*dlm3 + 126.28769397184625*dlm) + omx3*(-124.10987269781114*dlm2 + 27.429976953362896*dlm3 + 404.8743868729487*dlm) + (-12.362585066214681*dlm2 + 0.00013935927496978607*dlm3 + 70.98480615207401*dlm)*omx + (4889.000322137329 + 441.28447966969014*dl2 + 70.19038282912722*dl3 - 1798.0939409073276*dl)*x) + CA*d33c*_nf*(2.789630128993521e-6 + 4921.928178511819*x4 - 1937.2737574599053*x5 + 572.3341159882642*x6 - 126.34143068647165*x7 + 14.953465495510972*x8 + x2*(4112.389215159002 - 4794.148616796044*dl2 - 2184.6028365689886*dl3 - 437.4971457920455*dl4 - 16595.40050162117*dl) + x3*(4385.899382108362 - 2263.80240047488*dl2 + 815.6999149196427*dl3 + 112.45292920411022*dl4 - 4.024215883505673*dl) + omx3*(1054.6108655826722*dlm2 - 254.16191987627218*dlm3 - 4513.836078808586*dlm) + omx2*(-186.4834060522691*dlm2 - 10.551754754034768*dlm3 - 1393.3402344655901*dlm) + (-10.001790392395147*dlm2 - 0.0016132972331901125*dlm3 - 586.2047281592162*dlm)*omx + (-28346.73119239501 - 1694.0538250599097*dl2 - 705.1757403326384*dl3 + 78.8935012478646*dl4 - 1180.4610053078118*dl)*x) + CF*d33c*_nf*(-2.8938674895671293e-6 + 5725.068868393673*x4 - 1304.1365405152335*x5 + 353.590048501875*x6 - 82.21475727039513*x7 + 13.717098299273657*x8 + x3*(-8733.584011596065 + 13484.017506557846*dl2 - 3195.326081487795*dl3 + 488.77155721924305*dl4 - 31860.21164357397*dl) + x2*(-8692.735834158228 + 17170.373971256962*dl2 + 5433.814008202391*dl3 + 599.629735221736*dl4 + 35591.735712496644*dl) + omx2*(137.6122912642083*dlm2 + 5.31431824468334*dlm3 + 520.5775210297081*dlm) + omx3*(-1554.6381462366246*dlm2 + 163.07101012866482*dlm3 + 4795.342783328088*dlm) + (22.70277878442328*dlm2 + 0.0004253556673694253*dlm3 + 576.8211260294426*dlm)*omx + (26598.093426532727 + 2782.211975825979*dl2 + 868.7771646065672*dl3 - 157.84177341388818*dl4 + 3581.392736654987*dl)*x) + nf2*d33c*(-4522.666666666667 - (64*dl4)/3. - (64*dl5)/45. + (6976*zeta2)/3. + dl3*(-42.666666666666664 + (128*zeta2)/9.) + dl2*(-1194.6666666666667 + 768*zeta2 - 256*zeta3) - 3200*zeta3 + 1280*zeta2*zeta3 + dl*(-2944 - (1664*zeta2)/3. + (4096*zeta3)/3. - (320*zeta4)/3.) - (1664*zeta4)/3. - (4864*zeta5)/3.) + CA*d33c*_nf*(-97834.66666666667 + (256*dl5)/45. - (16*dl6)/5. - 6176*zeta32 + 44032*zeta2 + dl4*(-234.66666666666666 + 144*zeta2) + 52064*zeta3 - 15232*zeta2*zeta3 + dl3*(-1472 + (640*zeta2)/9. + 896*zeta3) + (26432*zeta4)/3. + dl2*(-9557.333333333334 + 3648*zeta2 + 4160*zeta3 + 1712*zeta4) + (28672*zeta5)/3. + dl*(-38506.666666666664 + 19968*zeta2 + (52352*zeta3)/3. - 7936*zeta2*zeta3 + 1760*zeta4 + 11136*zeta5) - 1160*zeta6) + CF*d33c*_nf*(18432 - (128*dl5)/15. + (32*dl6)/15. + 3840*zeta32 + dl4*(282.6666666666667 - (608*zeta2)/3.) - 36800*zeta2 + dl3*(725.3333333333334 + (512*zeta2)/3. - (2048*zeta3)/3.) - 13952*zeta3 + 5888*zeta2*zeta3 + 26288*zeta4 + dl2*(7936 - 8480*zeta2 - 512*zeta3 + 5024*zeta4) - 128*zeta5 + dl*(13312 - 15104*zeta2 - 12544*zeta3 + 7680*zeta2*zeta3 + 1920*zeta4 + 768*zeta5) - (1520*zeta6)/3.);
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
