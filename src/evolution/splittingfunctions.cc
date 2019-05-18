//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/splittingfunctions.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

namespace apfel
{
  /**
   * @brief The LO space-like splitting function classes
   */
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
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }

  //_________________________________________________________________________________
  P0qg::P0qg():
    Expression()
  {
  }
  double P0qg::Regular(double const& x) const
  {
    return 4 * TR * ( 1 - 2 * x + 2 * x * x );
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
    return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
  }

  /**
   * @brief The NLO space-like splitting function classes
   */
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
    const double ln1mx = log(1-x);
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
    const double x1nspa = gqq1 - gqq1l;
    return x1nspa;
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
    const double x1nsc = log(1-x) * _a2 + p1delta;
    return x1nsc;
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
    const double ln1mx = log(1-x);
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
    const double x1nsma = gqq1 - gqq1l;
    return x1nsma;
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
    const double x1psa =
      _nf * CF * ( - 8 + 24 * x - 224 / 9. * x * x + 80 / 9. / x
                   + ( 4 + 20 * x ) * lnx + 32 / 3. * x * x * lnx
                   - ( 4 + 4 * x ) * lnx2 );
    return x1psa;
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
    const double ln1mx = log(1-x);
    const double pqg   = x * x + ( 1 - x ) * ( 1 - x );
    const double pqgmx = x * x + ( 1 + x ) * ( 1 + x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double x1qga =
      + 2 * CF * _nf * ( 4  + 4 * ln1mx + ( 10 - 4 * ( ln1mx - lnx ) + 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 2 * Pi2 / 3 ) * pqg
                         - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) - 9 * x )
      + 2 * CA * _nf * ( 182 / 9. - 4 * ln1mx
                         + ( - 218 / 9. + 4 * ln1mx - 2 * ln1mx * ln1mx + 44 * lnx / 3 - lnx2 + Pi2 / 3 ) * pqg
                         + 2 * pqgmx * S2x + 40 / ( 9 * x ) + 14 * x / 9 - lnx2 * ( 2 + 8 * x )
                         + lnx * ( - 38 / 3. + 136 * x / 3 ) );
    return x1qga;
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
    const double ln1mx = log(1-x);
    const double pgq   = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
    const double pgqmx = - ( 1 + ( 1 + x ) * ( 1 + x ) ) / x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double x1gqa =
      + 2 * CF * _nf * ( - ( 20 / 9. + 4 * ln1mx / 3 ) * pgq - 4 * x / 3 )
      + 4 * CF * CF  * ( - 2.5 - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq - lnx2 * ( 1 - x / 2 ) - 7 * x / 2
                         - 2 * ln1mx * x + lnx * ( 2 + 7 * x / 2 ) )
      + 4 * CA * CF  * ( 28 / 9. + pgq * ( 0.5 + 11 * ln1mx / 3 + ln1mx * ln1mx - 2 * ln1mx * lnx + lnx2 / 2 - Pi2 / 6 ) + pgqmx * S2x
                         + 65 * x / 18 + 2 * ln1mx * x + 44 * x * x / 9 + lnx2 * ( 4 + x ) - lnx * ( 12 + 5 * x + 8 * x * x / 3 ) );
    return x1gqa;
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
    const double ln1mx = log(1-x);
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
    const double x1gga = ggg1 - ggg1l;
    return x1gga;
  }
  double P1gg::Singular(double const& x) const
  {
    return _a2g / ( 1 - x );
  }
  double P1gg::Local(double const& x) const
  {
    const double p1delta = ( - 2 * CF - 8 / 3. * CA ) * _nf + ( 32 / 3. + 12 * zeta3 ) * CA * CA;
    const double x1ggc = log(1-x) * _a2g + p1delta;
    return x1ggc;
  }

  /**
   * @brief The NNLO space-like splitting function classes (parametrized)
   */
  //_________________________________________________________________________________
  P2nsp::P2nsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2nsp::Regular(double const& x) const
  {
    const double x_2    = x * x;
    const double x_3    = x_2 * x;
    const double dl     = log(x);
    const double dl_2   = dl * dl;
    const double dl_3   = dl_2 * dl;
    const double dl_4   = dl_3 * dl;
    const double dl1    = log(1-x);
    const double d81    = 1. / 81.;
    const double p2nspa =
      1641.1 - 3135. * x + 243.6 * x_2 - 522.1 * x_3
      + 128. * d81 * dl_4 + 2400. * d81 * dl_3
      + 294.9 * dl_2 + 1258. * dl
      + 714.1 * dl1 + dl * dl1 * ( 563.9 + 256.8 * dl )
      + _nf * ( -197.0 + 381.1 * x + 72.94 * x_2 + 44.79 * x_3
                - 192. * d81 * dl_3  - 2608. * d81 * dl_2 - 152.6 * dl
                - 5120. * d81 * dl1 - 56.66 * dl * dl1 - 1.497 * x * dl_3 )
      + _nf * _nf * ( 32. * x * dl / ( 1 - x ) * ( 3. * dl + 10. ) + 64.
                      + ( 48. * dl_2 + 352. * dl + 384. ) * ( 1 - x ) ) * d81;
    return p2nspa;
  }
  double P2nsp::Singular(double const& x) const
  {
    const double p2nsb = ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1 - x );
    return p2nsb;
  }
  double P2nsp::Local(double const& x) const
  {
    const double dl1    = log(1-x);
    const double p2nspc =
      1174.898 * dl1 + 1295.624 - 0.24
      - _nf * ( 183.187 * dl1 + 173.938 - 0.011 )
      + _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
    return p2nspc;
  }

  //_________________________________________________________________________________
  P2nsm::P2nsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2nsm::Regular(double const& x) const
  {
    const double x_2    = x * x;
    const double x_3    = x_2 * x;
    const double dl     = log(x);
    const double dl_2   = dl * dl;
    const double dl_3   = dl_2 * dl;
    const double dl1    = log(1-x);
    const double d81    = 1. / 81.;
    const double p2nsma =
      1860.2 - 3505.* x + 297.0 * x_2 - 433.2 * x_3
      + 116. * d81 * dl_3 * dl + 2880. * d81 * dl_3
      + 399.2 * dl_2 + 1465.2 * dl
      + 714.1 * dl1 + dl * dl1 * ( 684.0 + 251.2 * dl )
      + _nf * ( -216.62 + 406.5 * x + 77.89 * x_2 + 34.76 * x_3
                - 256. * d81 * dl_3  - 3216. * d81 * dl_2 - 172.69 * dl
                - 5120. * d81 * dl1 - 65.43 * dl * dl1 - 1.136 * x * dl_3 )
      + _nf * _nf * ( 32.* x * dl / ( 1 - x ) * ( 3. * dl + 10. ) + 64.
                      + ( 48.* dl_2 + 352.* dl + 384. ) * ( 1.-x ) ) * d81;
    return p2nsma;
  }
  double P2nsm::Singular(double const& x) const
  {
    const double p2nsb = ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1 - x );
    return p2nsb;
  }
  double P2nsm::Local(double const& x) const
  {
    const double dl1    = log(1-x);
    const double p2nsmc =
      1174.898 * dl1 + 1295.624 - 0.154
      - _nf * ( 183.187 * dl1 + 173.938  - 0.005 )
      + _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
    return p2nsmc;
  }

  //_________________________________________________________________________________
  P2nss::P2nss(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2nss::Regular(double const& x) const
  {
    const double x_2    = x * x;
    const double d27    = 1. / 27.;
    const double dl     = log(x);
    const double dl_2   = dl * dl;
    const double dl_3   = dl_2 * dl;
    const double dl_4   = dl_3 * dl;
    const double x1     = 1 - x;
    const double dl1    = log(x1);
    const double p2nssa =
      x1 * ( 151.49 + 44.51 * x - 43.12 * x_2 + 4.820 * x_2 * x )
      + 40. * d27 * dl_4 - 80. * d27 * dl_3 + 6.892 * dl_2
      + 178.04 * dl + dl * dl1 * ( - 173.1 + 46.18 * dl )
      + x1 * dl1 * ( - 163.9 / x - 7.208 * x );
    return _nf * p2nssa;
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
    const double dl1   = log(1-x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double  p2ps1 =
      - 3584. / ( 27. * x ) * dl - 506.0 / x + 160. / 27. * dl_4
      - 400. / 9. * dl_3 + 131.4 * dl_2 - 661.6 * dl
      - 5.926  * dl1_3 - 9.751 * dl1_2 - 72.11 * dl1
      + 177.4 + 392.9 * x - 101.4 * x_2 - 57.04 * dl * dl1;
    const double p2ps2  =
      256. / ( 81. * x ) + 32. / 27. * dl_3 + 17.89 * dl_2
      + 61.75 * dl + 1.778 * dl1_2 + 5.944 * dl1 + 100.1
      - 125.2 * x + 49.26 * x_2 - 12.59 * x_3
      - 1.889 * dl * dl1;
    const double p2psa = ( 1 - x ) * _nf * ( p2ps1 + _nf * p2ps2 );
    return p2psa;
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
    const double dl1   = log(1-x);
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
    const double p2qga = _nf * ( p2qg1 + _nf * p2qg2 );
    return p2qga;
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
    const double dl1   = log(1-x);
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
    const double p2gqa = ( p2gq0 + _nf * ( p2gq1 + _nf * p2gq2 ) );
    return p2gqa;
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
    const double dl1   = log(1-x);
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
    const double p2gga = p2gga0 + _nf * ( p2gga1 + _nf * p2gga2 );
    return p2gga;
  }
  double P2gg::Singular(double const& x) const
  {
    const double p2ggb = ( 2643.521 - _nf * 412.172 - _nf * _nf * 16. / 9. ) / ( 1 - x );
    return p2ggb;
  }
  double P2gg::Local(double const& x) const
  {
    const double dl1   = log(1-x);
    const double p2ggc =
      2643.521 * dl1 + 4425.448 + 0.446
      - _nf * ( 412.172 * dl1 +  528.720 + 0.003 )
      + _nf * _nf * ( - 16. / 9. * dl1 + 6.4630 );
    return p2ggc;
  }

  /**
   * @brief The NNNLO splitting function classes (parametrized and
   * leading color). Only the +, -, and valence contributions have
   * been computed so far.
   *
   */
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

    // Leading large-n_c, nf^0 and nf^1, parametrized.
    const double p3nsa0 =
      2.5e4 * ( omy * ( 3.5254 + 8.6935 * y - 1.5051 * y2 + 1.8300 * y3 )
                + 11.883 * y * dl - 0.09066 * y * dl2 + 11.410 * omy * dlm + 13.376 * dl * dlm )
      + 5.167133e4 * dl + 1.712095e4 * dl2 + 2.863226e3 * dl3 + 2.978255e2 * dl4
      + 1.6e1 * dl5 + 5e-1 * dl6 - 2.973385e4 + 1.906980e4 * dlm;
    const double p3nsa1 =
      2.5e4* ( omy * ( - 0.74077 + 1.4860 * y - 0.23631 * y2 + 0.31584 * y3 )
               + 2.5251 * omy * dlm + 2.5203 * dl * dlm + 2.2242 * y * dl
               - 0.02460 * y * dl2 + 0.00310 * y * dl3 )
      - 9.239374e3 * dl - 2.917312e3 * dl2 - 4.305308e2 * dl3 - 3.6e1 * dl4
      - 4. / 3. * dl5 + 8.115605e3 - 3.079761e3 * dlm;

    // Nonleading large-n_c, nf^0 and nf^1: two approximations
    const double p3npa01 =
      3948.16 * omy - 2464.61 * ( 2 * y - y2 ) * omy - 1839.44 * dl2 - 402.156 * dl3
      - 1777.27 * dlm2 * omy - 204.183 * dlm3 * omy + 507.152 - 5.587553e+1 * dl4
      - 2.831276 * dl5 - 1.488340e-1 * dl6 - 2.601749e+3 - 2.118867e+3 * dlm;
    const double p3npa02 =
      ( 8698.39 - 10490.47 * y ) * y * omy + 1389.73 * dl + 189.576 * dl2
      - 173.936 * dlm2 * omy + 223.078 * dlm3 * omy + 505.209 - 5.587553e+1 * dl4
      - 2.831276 * dl5 - 1.488340e-1 * dl6 - 2.601749e+3 - 2.118867e+3 * dlm;

    const double p3npa11 =
      ( - 1116.34 + 1071.24 * y ) * y * omy - 59.3041 * dl2 - 8.4620 * dl3
      - 143.813 * dlm * omy - 18.8803 * dlm3 * omy - 7.33927 + 4.658436 * dl4
      + 2.798354e-1 * dl5 + 3.121643e+2 + 3.379310e+2 * dlm;
    const double p3npa12 =
      ( - 690.151 - 656.386 * y2 ) * omy + 133.702 * dl2 + 34.0569 * dl3
      - 745.573 * dlm * omy + 8.61438 * dlm3 * omy - 7.53662 + 4.658437 * dl4
      + 2.798354e-1 * dl5 + 3.121643e+2 + 3.379310e+2 * dlm;

    // nf^2 (parametrized) and nf^3 (exact)
    const double p3nspa2 =
      2.5e2 *  ( omy * ( 3.0008 + 0.8619 * y - 0.12411 * y2 + 0.31595 * y3 )
                 - 0.37529 * y * dl - 0.21684 * y * dl2 - 0.02295 * y * dl3
                 + 0.03394 * omy * dlm + 0.40431  * dl * dlm )
      + 3.930056e+2 * dl + 1.125705e+2 * dl2 + 1.652675e+1 * dl3
      + 7.901235e-1 * dl4 - 3.760092e+2 + 2.668861e+1 * dlm;
    const double p3nsa3  =
      - 2.426296e0 - 8.460488e-1 * y + ( 5.267490e-1 * dm - 3.687243e0 + 3.160494e0 * y ) * dl
      - ( 1.316872e0 * ( dm + 1e-1 ) - 1.448560e0 * y ) * dl2
      - ( 2.633745e-1 * dm - 1.31687e-1 * ( 1 + y ) ) * dl3;

    // Assembly
    const double p3nspai = p3nsa0 + _nf * ( p3nsa1 + _nf * ( p3nspa2 + _nf * p3nsa3 ) );
    if (_imod == 1)
      return p3nspai + p3npa01 + _nf * p3npa11;
    else if (_imod == 2)
      return p3nspai + p3npa02 + _nf * p3npa12;
    else
      return p3nspai + 0.5* ( ( p3npa01 + p3npa02 ) + _nf * ( p3npa11 + p3npa12 ) );
  }
  double P3nsp::Singular(double const& y) const
  {
    const double d1 = 1 / ( 1 - y );
    const double a4qi =
      2.120902e+4
      - 5.179372e+3 * _nf
      + 1.955772e+2 * _nf * _nf
      + 3.272344 * _nf * _nf * _nf;
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
    const double dl1 = log(1-y);
    const double a4qi  =
      2.120902e+4
      - 5.179372e+3* _nf
      + 1.955772e+2* _nf * _nf
      + 3.272344* _nf * _nf * _nf;
    const double a4ap1 = - 507.152 + 7.33927 * _nf;
    const double a4ap2 = - 505.209 + 7.53662 * _nf;

    const double b4qi =
      2.579609e+4 + 0.08
      - ( 5.818637e+3 + 0.97 ) * _nf
      + ( 1.938554e+2 + 0.0037)* _nf * _nf
      + 3.014982 * _nf * _nf * _nf;
    const double b4ap1 = - 2405.03 + 267.965 * _nf;
    const double b4ap2 = - 2394.47 + 269.028 * _nf;

    if (_imod == 1)
      return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap1;
    else if (_imod == 2)
      return ( a4qi + a4ap1 ) * dl1 + b4qi + b4ap2;
    else
      return ( a4qi + 0.5 * ( a4ap1 + a4ap2 ) ) * dl1
             + b4qi + 0.5 * ( b4ap1 + b4ap2 );
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

    const double p3nsa0  =
      2.5e+4 * ( omy * ( 3.5254 + 8.6935 * y - 1.5051 * y2 + 1.8300 * y3 )
                 + 11.883 * y * dl - 0.09066 * y * dl2 + 11.410 * omy * dlm + 13.376  * dl * dlm )
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
      ( 4043.59 - 15386.6 * y ) * y * omy + 502.481 + 1532.96  * dl2 + 31.6023 * dl3
      - 3997.39  * dlm * omy + 511.567 * dlm3 * omy + 4.964335e-1 * ( dl6 + 18 * dl5 )
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
                 + 2.7511 * omy * dlm + 3.2709  * dl * dlm )
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
    const double dl1 = log(1-y);

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
      + ( 1.938554e+2 + 0.0037)* _nf * _nf
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
      + 533.44 * dl4 + 97.27 *  dl5 + 4 * dl6 + 60.40 * omy * dlm2 + 4.685 * omy * dlm3;
    const double p3nsa12 =
      1030.79 * omy * y + 1266.77 * omy * ( 2 - y2 ) + 2987.83 * dl + 273.05 * dl2 - 923.48 * dl3
      - 236.76 * dl4 - 33.886 * dl5 - 4 * dl6 - 254.63 * omy * dlm - 0.28953 * omy * dlm3;

    // nf^2 (parametrized)
    const double p3nssa2 =
      2.5e+2 * ( omy * ( - 4.7656 + 1.6908 * y + 0.1703 * y2 )
                 - 0.41652 * y *dl + 0.90777 * y * dl2 + 0.12478 * y * dl3
                 + 0.17155 * omy * dlm + 0.17191  * dl * dlm )
      - 6.473971e+2 * dl - 6.641219e+1 * dl2 - 5.353347 * dl3 - 5.925926 *dl4
      - 3.950617e-1 * dl5 + 1.970002e+1 * omy * dlm - 3.435474 * omy * dlm2;

    if (_imod == 1)
      return _nf * p3nsa11 + _nf * _nf * p3nssa2;
    else if (_imod == 2)
      return _nf * p3nsa12 + _nf * _nf * p3nssa2;
    else
      return 0.5 *_nf * ( p3nsa11 + p3nsa12 ) + _nf * _nf * p3nssa2;
  }

  /**
   * @brief The LO space-like splitting function for tranversely
   * polarised PDFs. Reference
   * https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
  //_________________________________________________________________________________
  P0transns::P0transns():
    Expression()
  {
  }
  double P0transns::Regular(double const&) const
  {
    return - 4 * CF;
  }
  double P0transns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double P0transns::Local(double const& x) const
  {
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }

  /**
   * @brief The NLO space-like splitting function for tranversely
   * polarised PDFs. Reference
   * https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
  //_________________________________________________________________________________
  P1transnsp::P1transnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = 2 * CA * CF * ( 134. / 9. - 2 * Pi2 / 3 ) - 80 * _nf * CF * TR / 9;
  }
  double P1transnsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double func  = 2 * x * lnx / ( 1 - x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double gqq1  =
      + 4 * CF * CF * ( 1 - x - ( 3. / 2. + 2 * ln1mx ) * func )
      + 2 * CF * CA * ( - 143. / 9. + 2 * Pi2 / 3 + x + ( 11. / 3. + lnx ) * func )
      + 8 * _nf * CF * TR * ( - func + 10. / 3. ) / 3
      + 4 * CF * ( CF - CA / 2 ) * ( - 1 + x - 4 * S2x / ( 1 + x ) );
    return gqq1;
  }
  double P1transnsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1transnsp::Local(double const& x) const
  {
    const double p1delta =
      + 4 * CF * CF * ( 3. / 8. - Pi2 / 2 + 6 * zeta3 )
      + 2 * CF * CA * ( 17. / 12. + 11 * Pi2 / 9 - 6 * zeta3 )
      - 8 * _nf * CF * TR * ( 1. / 4. + Pi2 / 3 );
    const double x1nsc = log(1-x) * _a2 + p1delta;
    return x1nsc;
  }

  //_________________________________________________________________________________
  P1transnsm::P1transnsm(int const& nf):
    P1transnsp(nf)
  {
  }
  double P1transnsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double func  = 2 * x * lnx / ( 1 - x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double gqq1  =
      + 4 * CF * CF * ( 1 - x - ( 3. / 2. + 2 * ln1mx ) * func )
      + 2 * CF * CA * ( - 143. / 9. + 2 * Pi2 / 3 + x + ( 11. / 3. + lnx ) * func )
      + 8 * _nf * CF * TR * ( - func + 10. / 3. ) / 3
      - 4 * CF * ( CF - CA / 2 ) * ( - 1 + x - 4 * S2x / ( 1 + x ) );
    return gqq1;
  }
}
