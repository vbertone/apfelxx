//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionspol_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

namespace apfel
{
  /**
   * @brief The LO space-like longitudinally polarised splitting
   * function classes
   */
  //_________________________________________________________________________________
  P0polns::P0polns():
    P0ns()
  {
  }

  //_________________________________________________________________________________
  P0polqg::P0polqg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0polqg::Regular(double const& x) const
  {
    return 4 * _nf * TR * ( 2 * x - 1 );
  }

  //_________________________________________________________________________________
  P0polgq::P0polgq():
    Expression()
  {
  }
  double P0polgq::Regular(double const& x) const
  {
    return 2 * CF * ( 2 - x );
  }

  //_________________________________________________________________________________
  P0polgg::P0polgg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0polgg::Regular(double const& x) const
  {
    return 4 * CA * ( - 2 * x + 1 );
  }
  double P0polgg::Singular(double const& x) const
  {
    return 4 * CA / ( 1 - x );
  }
  double P0polgg::Local(double const& x) const
  {
    return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
  }

  /**
   * @brief The NLO space-like longitudinally polarised splitting
   * function classes
   */
  //_________________________________________________________________________________
  P1polnsp::P1polnsp(int const& nf):
    P1nsm(nf)
  {
  }

  //_________________________________________________________________________________
  P1polnsm::P1polnsm(int const& nf):
    P1nsp(nf)
  {
  }

  //_________________________________________________________________________________
  P1polps::P1polps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1polps::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return 8 * _nf * CF * TR * ( ( 1 - x ) - ( 1 - 3 * x ) * lnx - ( 1 + x ) * lnx2 );
  }

  //_________________________________________________________________________________
  P1polqg::P1polqg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1polqg::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double dpqg   = 2 * x - 1;
    const double dpqgmx = - 2 * x - 1;
    const double S2x    = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * _nf * TR * CF * ( - 22 + 27 * x - 9 * lnx + 8 * ( 1 - x ) * ln1mx + dpqg * ( 2 * ln1mx2 - 4 * ln1mx * lnx + lnx2 - 4 * zeta2 ) )
      + 4 * _nf * TR * CA * ( ( 24 - 22 * x ) - 8 * ( 1 - x ) * ln1mx + ( 2 + 16 * x ) * lnx - 2 * ( ln1mx2 - zeta2 ) * dpqg - ( 2 * S2x - 3 * lnx2 ) * dpqgmx );
  }

  //_________________________________________________________________________________
  P1polgq::P1polgq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1polgq::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double dpgq   = 2 - x;
    const double dpgqmx = 2 + x;
    const double S2x    = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * _nf * CF * TR * ( - 4 * ( x + 4 ) / 9. - 4 * dpgq * ln1mx / 3. )
      + 4 * CF * CF * ( - 1 / 2. - ( 4 - x ) * lnx / 2 - dpgqmx * ln1mx + ( - 4 - ln1mx2 + lnx2 / 2 ) * dpgq )
      + 4 * CF * CA * ( ( 4 - 13 * x ) * lnx + ( 10 + x ) * ln1mx / 3 + ( 41 + 35 * x ) / 9
                        + ( - 2 * S2x + 3 * lnx2 ) * dpgqmx / 2 + ( ln1mx2 - 2 * ln1mx * lnx - zeta2 ) * dpgq );
  }

  //_________________________________________________________________________________
  P1polgg::P1polgg(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2g = - 80 / 9. * CA * TR * nf + ( 268 / 9. - 8 * zeta2 ) * CA * CA;
  }
  double P1polgg::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1 - x);
    const double dpgg   = 1 / ( 1 - x ) - 2 * x + 1;
    const double dpggmx = 1 / ( 1 + x ) + 2 * x + 1;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double ggg1  =
      - 4 * CA * TR * _nf * ( 4 * ( 1 - x ) + 4 * ( 1 + x ) / 3. * lnx + 20 / 9. * dpgg )
      - 4 * CF * TR * _nf * ( 10 * ( 1 - x ) + 2 * ( 5 - x ) * lnx + 2 * ( 1 + x ) * lnx2 )
      + 4 * CA * CA * ( ( 29 - 67 * x ) * lnx / 3. - 19 * ( 1 - x ) / 2. + 4 * ( 1 + x ) * lnx2
                        - 2 * S2x * dpggmx + ( 67 / 9. - 4 * ln1mx * lnx + lnx2 - 2 * zeta2 ) * dpgg );
    const double ggg1l = _a2g / ( 1 - x );
    return ggg1 - ggg1l;
  }
  double P1polgg::Singular(double const& x) const
  {
    return _a2g / ( 1 - x );
  }
  double P1polgg::Local(double const& x) const
  {
    const double p1delta = ( - 2 * CF - 8 / 3. * CA ) * _nf + ( 32 / 3. + 12 * zeta3 ) * CA * CA;
    return log(1-x) * _a2g + p1delta;
  }

  /**
   * @brief The NNLO space-like longitudinally polarised splitting
   * function classes (parametrized)
   */
  //_________________________________________________________________________________
  P2polnsp::P2polnsp(int const& nf):
    P2nsm(nf)
  {
  }

  //_________________________________________________________________________________
  P2polnsm::P2polnsm(int const& nf):
    P2nsp(nf)
  {
  }

  //_________________________________________________________________________________
  P2polnss::P2polnss():
    Expression()
  {
  }

  //_________________________________________________________________________________
  P2polps::P2polps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2polps::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double x3     = x * x2;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    const double lnx4   = lnx * lnx3;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double P2ps1 = - 344./27. * lnx4 - (90.9198 + 81.50* x)* lnx3 - (368.6 - 349.9* x)* lnx2 - (739.0 - 232.57* ln1mx)* lnx
                         - 1362.6 + 1617.4 * x - 674.8 * x2 + 167.41 * x3 - 204.76 * ln1mx - 12.61 * ln1mx2 - 6.541 * ln1mx3;
    const double P2ps2 = (1.1741 - 0.8253* x)* lnx3  + (13.287 + 10.657* x)* lnx2 + 45.482 * lnx + 49.13 - 30.77 * x - 4.307 * x2
                         - 0.5094 *x3 + 9.517 * ln1mx + 1.7805 * ln1mx2;
    return ( 1 - x ) * _nf * ( P2ps1 + _nf * P2ps2 );
  }

  //_________________________________________________________________________________
  P2polqg::P2polqg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2polqg::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double x3     = x * x2;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    const double lnx4   = lnx * lnx3;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double ln1mx4 = ln1mx * ln1mx3;
    const double P2qg1 = - 151./3. * lnx4 - (385.64 + 73.30* x)* lnx3 - (894.8 - 1145.3* x)* lnx2 - (1461.2 - 825.4* ln1mx)* lnx
                         - 2972.4 + 4672.* x - 1221.6 * x2 - 18.0 * x3 + 278.32* ln1mx - 90.26* ln1mx2 - 5.30* ln1mx3 + 3.784*ln1mx4;
    const double P2qg2 = 16./9. * lnx4 + (30.739  + 10.186* x) * lnx3 + (196.96 + 179.1* x)* lnx2 + (526.3  - 47.30* ln1mx)* lnx
                         + 499.65 - 432.18 * x - 141.63 * x2 - 11.34 * x3 - 6.256 * ln1mx + 7.32 * ln1mx2 + 0.7374 * ln1mx3;
    return _nf * ( P2qg1 + _nf * P2qg2 );
  }

  //_________________________________________________________________________________
  P2polgq::P2polgq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2polgq::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double x3     = x * x2;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    const double lnx4   = lnx * lnx3;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double ln1mx4 = ln1mx * ln1mx3;
    const double P2gq0 = 11512./81.* lnx4 + (888.003  + 175.1* x)* lnx3 + (2140. - 850.7* x)* lnx2 + (4046.6 - 1424.8* ln1mx)* lnx
                         + 6159. - 3825.9 * x + 1942.* x2 - 742.1 * x3 + 1843.7* ln1mx + 451.55* ln1mx2 + 59.3* ln1mx3 + 5.143* ln1mx4;
    const double P2gq1 = - 128./27. * lnx4 - (39.3872 + 30.023*x)* lnx3 - (202.46 + 126.53* x)* lnx2 - (308.98 + 16.18* ln1mx)* lnx
                         - 301.07 - 296.0 * x + 406.13 * x2 - 101.62 * x3 - 171.78* ln1mx - 47.86 * ln1mx2 - 4.963 * ln1mx3;
    const double P2gq2 = 16./27. * ( - 12. + 10.* x + ( 8.+ 2.*x)* ln1mx + (6.- 3.*x)* ln1mx2 );
    return P2gq0 + _nf * ( P2gq1 + _nf * P2gq2 );
  }

  //_________________________________________________________________________________
  P2polgg::P2polgg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2polgg::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double x3     = x * x2;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    const double lnx4   = lnx * lnx3;
    const double ln1mx  = log(1 - x);
    const double P2ggA0 = 504. * lnx4 + (3777.5  + 1167.* x)* lnx3 + (10902. - 863.* x)* lnx2 + (23091. - 12292.* ln1mx)* lnx
                          + 30988. - 39925.* x + 13447.* x2 - 4576.* x3 - 13247.* (1.-x)*ln1mx + 3801.* ln1mx;
    const double P2ggA1 = - 766./27. * lnx4 - (357.798 - 131.* x)* lnx3 - (1877.2 - 613.1* x)* lnx2 - (3524. + 7932.* ln1mx)* lnx
                          - 1173.5 + 2648.6 * x - 2160.8 * x2 + 1251.7 * x3 - 6746.* (1.-x)*ln1mx - 295.7* ln1mx;
    const double P2ggA2 = - 1.1809 * lnx3 - (6.679 - 15.764* x)* lnx2 - (13.29 + 16.944* ln1mx) * lnx - 16.606 + 32.905 * x
                          - 18.30 * x2 + 2.637 * x3 - 0.210 * ln1mx;
    return P2ggA0 + _nf * ( P2ggA1 + _nf * P2ggA2 );
  }
  double P2polgg::Singular(double const& x) const
  {
    return ( 2643.521 - _nf * ( 412.172 + _nf * 16. / 9. ) ) / ( 1 - x );
  }
  double P2polgg::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return 2643.521 * ln1mx + 4425.448 + 2.314 - _nf * ( 412.172 * ln1mx + 528.720 - 0.184 ) - _nf * _nf * ( 16. / 9. * ln1mx - 6.4630 + 0.0023 );
  }
}
