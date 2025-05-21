//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionsunp_sl_qed.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  P01qedns::P01qedns():
    Expression()
  {
  }
  double P01qedns::Regular(double const& x) const
  {
    return - 2 * ( 1 + x );
  }
  double P01qedns::Singular(double const& x) const
  {
    return 4 / ( 1 - x );
  }
  double P01qedns::Local(double const& x) const
  {
    return 4 * log(1 - x) + 3;
  }

  //_________________________________________________________________________________
  P01qedqgm::P01qedqgm():
    Expression()
  {
  }
  double P01qedqgm::Regular(double const& x) const
  {
    return 4 * ( 1 - 2 * x + 2 * x * x );
  }

  //_________________________________________________________________________________
  P01qedgmq::P01qedgmq():
    Expression()
  {
  }
  double P01qedgmq::Regular(double const& x) const
  {
    return 4 * ( - 1 + 0.5 * x + 1 / x );
  }

  //_________________________________________________________________________________
  P01qedgmgm::P01qedgmgm():
    Expression()
  {
  }
  double P01qedgmgm::Local(double const&) const
  {
    return - 4. / 3;
  }

  //_________________________________________________________________________________
  P02qednsp::P02qednsp(double const& crat2):
    Expression(),
    _crat2(crat2)
  {
    _a2 = - 80 * _crat2 / 9.;
  }
  double P02qednsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1 + x) - Pi2 / 6;
    const double gqq1  =
      + 4 * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x ) - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
      + 4 * _crat2 * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 )
      + 4 * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    return gqq1 - gqq1l;
  }
  double P02qednsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P02qednsp::Local(double const& x) const
  {
    const double p1delta = 3 / 2. + 24 * zeta3 - 2 * Pi2 + _crat2 * ( - 2 / 3. - 8 * Pi2 / 9. );
    return log(1 - x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P02qednsm::P02qednsm(double const& crat2):
    Expression(),
    _crat2(crat2)
  {
    _a2 = - 80 * _crat2 / 9.;
  }
  double P02qednsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1 + x) - Pi2 / 6;
    const double gqq1  =
      + 4 * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x ) - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
      + 4 * _crat2 * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 )
      - 4 * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    return gqq1 - gqq1l;
  }
  double P02qednsm::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P02qednsm::Local(double const& x) const
  {
    const double p1delta = 3 / 2. + 24 * zeta3 - 2 * Pi2 + _crat2 * ( - 2 / 3. - 8 * Pi2 / 9. );
    return log(1 - x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P02qedps::P02qedps():
    Expression()
  {
  }
  double P02qedps::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return 2 * ( - 8 + 24 * x - 224 / 9. * x * x + 80 / 9. / x
                 + ( 4 + 20 * x ) * lnx + 32 / 3. * x * x * lnx
                 - ( 4 + 4 * x ) * lnx2 );
  }

  //_________________________________________________________________________________
  P02qedqgm::P02qedqgm():
    Expression()
  {
  }
  double P02qedqgm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqg   = x * x + ( 1 - x ) * ( 1 - x );
    return 4 * ( 4  + 4 * ln1mx + ( 10 - 4 * ( ln1mx - lnx ) + 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 2 * Pi2 / 3 ) * pqg
                 - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) - 9 * x );
  }

  //_________________________________________________________________________________
  P02qedgmq::P02qedgmq(double const& crat2):
    Expression(),
    _crat2(crat2)
  {
  }
  double P02qedgmq::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pgq   = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
    return 4 * ( - 2.5 - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq - lnx2 * ( 1 - x / 2 ) - 7 * x / 2
                 - 2 * ln1mx * x + lnx * ( 2 + 7 * x / 2 ) )
           + 4 * _crat2 * ( - ( 20 / 9. + 4 * ln1mx / 3 ) * pgq - 4 * x / 3 );
  }

  //_________________________________________________________________________________
  P02qedgmgm::P02qedgmgm():
    Expression()
  {
  }
  double P02qedgmgm::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return 4 * ( - 16 + 4 / ( 3 * x ) + 8 * x + ( 20 * x * x ) / 3 - lnx2 * ( 2 + 2 * x ) - lnx * ( 6 + 10 * x ) );
  }
  double P02qedgmgm::Local(double const&) const
  {
    return - 4.;
  }

  //_________________________________________________________________________________
  P11qednsp::P11qednsp():
    Expression()
  {
  }
  double P11qednsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1 + x) - Pi2 / 6;
    return 8 * CF * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - lnx * ( 1.5 + 7 * x / 2 ) - lnx2 * ( 1 + x ) / 2 - 5 * ( 1 - x ) )
           + 8 * CF *( 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) + 2 * pqqmx * S2x );
  }
  double P11qednsp::Local(double const&) const
  {
    return 8 * CF * ( 3 / 8. + 6 * zeta3 - Pi2 / 2 );
  }

  //_________________________________________________________________________________
  P11qednsm::P11qednsm():
    Expression()
  {
  }
  double P11qednsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1 + x) - Pi2 / 6;
    return 8 * CF * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - lnx * ( 1.5 + 7 * x / 2 ) - lnx2 * ( 1 + x ) / 2 - 5 * ( 1 - x ) )
           - 8 * CF *( 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) + 2 * pqqmx * S2x );
  }
  double P11qednsm::Local(double const&) const
  {
    return 8 * CF * ( 3 / 8. + 6 * zeta3 - Pi2 / 2 );
  }

  //_________________________________________________________________________________
  P11qedqg::P11qedqg():
    Expression()
  {
  }
  double P11qedqg::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqg   = x * x + ( 1 - x ) * ( 1 - x );
    return 4 * TR * ( 4 - 9 * x - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) + 4 * ln1mx + ( 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 4 * ( ln1mx - lnx ) - 2 * Pi2 / 3 + 10 ) * pqg );
  }

  //_________________________________________________________________________________
  P11qedqgm::P11qedqgm():
    Expression()
  {
  }
  double P11qedqgm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqg   = x * x + ( 1 - x ) * ( 1 - x );
    return 4 * CF * ( 4 - 9 * x - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) + 4 * ln1mx + ( 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 4 * ( ln1mx - lnx ) - 2 * Pi2 / 3 + 10 ) * pqg );
  }

  //_________________________________________________________________________________
  P11qedgq::P11qedgq():
    Expression()
  {
  }
  double P11qedgq::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pgq   = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
    return 4 * CF * ( - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq  + lnx * ( 2 + 7 * x / 2 ) - lnx2 * ( 1 - x / 2 ) - 2 * ln1mx * x  - 7 * x / 2 - 2.5 );
  }

  //_________________________________________________________________________________
  P11qedgmq::P11qedgmq():
    P11qedgq()
  {
  }

  //_________________________________________________________________________________
  P11qedgg::P11qedgg():
    Expression()
  {
  }
  double P11qedgg::Local(double const&) const
  {
    return - 4 * TR;
  }

  //_________________________________________________________________________________
  P11qedggm::P11qedggm():
    Expression()
  {
  }
  double P11qedggm::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return 4 * CF * ( - 16 + 8 * x + 20. * x2 / 3 + 4 / 3. / x - ( 6 + 10 * x ) * lnx - 2 * ( 1 +  x ) * lnx2 );
  }

  //_________________________________________________________________________________
  P11qedgmg::P11qedgmg():
    Expression()
  {
  }
  double P11qedgmg::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return 4 * TR * ( - 16 + 8 * x + 20. * x2 / 3 + 4 / 3. / x - ( 6 + 10 * x ) * lnx - 2 * ( 1 +  x ) * lnx2 );
  }

  //_________________________________________________________________________________
  P11qedgmgm::P11qedgmgm():
    Expression()
  {
  }
  double P11qedgmgm::Local(double const&) const
  {
    return - 4 * CF;
  }

  //_________________________________________________________________________________
  P21qednsp::P21qednsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P21qednsp::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double lx4  = lx * lx3;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x * l1x2;
    const double l1x4 = l1x * l1x3;
    const double kns20 = 128. / 27.* l1x4 + 112. / 9.* l1x3 + 175.3 * l1x2 + 142.3 * l1x + 1353. - 1262.* x + 449.2 * x2
                         - 1445.* x3 - lx * l1x * ( 162.7 * lx + 195.4 * l1x ) + 1169. * x * lx + 50.08 * ( 1 - x ) * l1x3 + 744.6 * lx
                         + 201.6 * lx2 + 80. / 3. * lx3 + 64. / 27.* lx4;
    const double kns21 = - 32. / 27.* l1x3 - 11.858 * l1x2 - 18.77 * l1x - 40.035 + 114.4 * x - 24.86 * x2 - 53.39 * x3
                         + lx * l1x * ( 8.523 * lx + 269.4 * l1x ) - 26.63 * x * lx + 270. * ( 1 - x ) * l1x2 - 21.55 * lx - 10.992 * lx2 - 32. / 27.* lx3;
    return kns20 + _nf * kns21;
  }

  //_________________________________________________________________________________
  P21qedps::P21qedps():
    Expression()
  {
  }
  double P21qedps::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double lx  = log(x);
    const double lx2 = lx * lx;
    const double lx3 = lx * lx2;
    const double lx4 = lx * lx3;
    return ( 2464./81./ x -  432. - 72.* x + 38360./81.* x2 - lx * ( 344. + 368.* x + 3584./27.* x2 )
             - lx2 * ( 144. + 104.* x + 224./9.* x2 ) - lx3 * ( 16. - 16.* x - 128./9.* x2 ) - lx4 * 8./3.* (1. - 2.* x) ) * 4./3.;
  }

  //_________________________________________________________________________________
  P21qedggm::P21qedggm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P21qedggm::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x * l1x2;
    const double kg20 = 32. / 27. * l1x3 - 79.13 * l1x2 + 87.22 * l1x + 1738. - 1580.* x - 160.* x2 - 566.7 * x3
                        - lx * l1x * ( 549.5 + 1230. * lx + 433.2 * l1x ) + 2176.* lx + 1123.7 * lx2 + ( 2400. + 448.* lx ) / 27. * lx3
                        - (73.1409 - 128./3. * lx) / x;
    const double kg21 = - 32. / 9.* l1x2 + 16.38 * l1x + 68.10 - 36.42 * x + 56.95 * x2 - 44.10 * x3
                        - lx* l1x * ( 16.18 + 38.33 * lx + 9.133 * l1x) - 10.76 * lx + 26.41 * lx2 - 64./27.* lx3 - 40.5597 / x;
    return ( 1 - x ) * ( kg20 + _nf * kg21 );
  }

  //_________________________________________________________________________________
  DeltaP11qednsp::DeltaP11qednsp():
    Expression()
  {
  }
  double DeltaP11qednsp::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnomx  = log(1 - x);
    const double lnomx2 = lnomx * lnomx;
    return 4. * CF * ( 7. - 10. * x - Pi2 / 6. * ( 6. - 12. * x + 16. * x2 ) + ( 1. - 16. * x + 32. * x2 ) * lnx
                       + ( 1. - 2. * x + 4. * x2 ) * lnx2 - ( 5. - 36. * x + 32. * x2 ) * lnomx
                       + ( 4. - 8 * x + 8. * x2 ) * ( lnomx2 - lnx * lnomx ) + ( 2. - 4. * x + 8 * x2 ) * dilog(x) );
  }

  //_________________________________________________________________________________
  DeltaP11qedggm::DeltaP11qedggm():
    Expression()
  {
  }
  double DeltaP11qedggm::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double lnomx = log(1 - x);
    return 8. * CF * ( 2. / 3. / x - 20. / 3. + 2. * x / 3. + 16. * x2 / 3. - ( 1. + 5. * x - 4. * x2 / 3. ) * lnx
                       - ( 1. + x ) * lnx2 + ( 4. / 3. / x + 1. - x - 4. * x2 / 3. ) * lnomx - 2. * ( 1. + x ) * ( dilog(x) - Pi2 / 6. ) );
  }

  //_________________________________________________________________________________
  DeltaP21qednsp::DeltaP21qednsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double DeltaP21qednsp::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double L1  = log(1 - x);
    const double L12 = L1 * L1;
    const double L13 = L1 * L12;
    const double L14 = L1 * L13;
    const double L0  = log(x);
    const double L02 = L0 * L0;
    const double L03 = L0 * L02;
    const double L04 = L0 * L03;
    return 9.482 * L14 + ( 33.37 + 2585 * ( 1. - x ) ) * L13 + 122.6 * L12 + ( 5598. + 7949. * ( 1. - x ) ) * L1
           - 9825. * L0 * L1 - 2.963 * L04 + 7.407 * L03 - ( 176.0 - 2616 * x ) * L02 - 828.6 * L0 - 1851. + 30120. * x - 7595. * x2
           + _nf * ( - ( 1.044 + 32.57 * ( 1. - x ) ) * L13 + 26.75 * L12 - ( 0.266 - 615.1 * ( 1. - x ) ) * L1 + 4.557 * L0 * L1
                     - 0.529 * L03 + ( 12.23 - 38.59 * x ) * L02 + 41.91 * L0 + 75.74 + 733.7 * x- 1003. * x2 );
  }
  double DeltaP21qednsp::Local(double const&) const
  {
    return - 0.65 + _nf * 0.05;
  }

  //_________________________________________________________________________________
  DeltaP21qedps::DeltaP21qedps():
    Expression()
  {
  }
  double DeltaP21qedps::Regular(double const& x) const
  {
    const double L1  = log(1 - x);
    const double L12 = L1 * L1;
    const double L0  = log(x);
    const double L02 = L0 * L0;
    const double L03 = L0 * L02;
    const double L04 = L0 * L03;
    return ( 2.083 * L12 - 68.96 * L1 ) * pow(1. - x, 3) - 86.40 * pow(1. - x, 2) * L0 * L1 + 1.778 * L04 + 3.278 * L03
           + ( 41.86 + 105.6 * x ) * ( 1. - x ) * L02 + 1.241 * pow(1. - x, 2) * L0- ( 15.80 / x + 94.39 - 8.281 * x ) * pow(1. - x, 3);
  }

  //_________________________________________________________________________________
  DeltaP21qedggm::DeltaP21qedggm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double DeltaP21qedggm::Regular(double const& x) const
  {
    const double L1  = log(1 - x);
    const double L12 = L1 * L1;
    const double L13 = L1 * L12;
    const double L0  = log(x);
    const double L02 = L0 * L0;
    const double L03 = L0 * L02;
    const double L04 = L0 * L03;
    return -( 25.77 * L13 + 766.7 * L1 ) * ( 1. - x ) - 1337. * L0 * L1 - 4.741 * L04 + 4.741 * L03
           + ( 83.11 - 464.9 * x ) * L02 + 443.1 * L0 - ( 60.36 / x - 1772. - 650.3 * x ) * ( 1. - x )
           + _nf * ( - ( 0.737 * L13 - 216. * L1 ) * ( 1. - x ) + 310.6 * L0 * L1 - ( 40.99 - 36.66 * x ) * L02
                     - 113.7 * L0 + ( 15.80 / x - 331.3 - 48.88 * x ) * ( 1. - x ) );
  }
}
