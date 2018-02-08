//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/splittingfunctions_tl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  /**
   * @brief The LO time-like splitting function classes
   */
  //_________________________________________________________________________________
  P0Tns::P0Tns():
    Expression()
  {
  }
  double P0Tns::Regular(double const& x) const
  {
    return - 2 * CF * ( 1 + x );
  }
  double P0Tns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double P0Tns::Local(double const& x) const
  {
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }
  
  //_________________________________________________________________________________
  P0Tqg::P0Tqg():
    Expression()
  {
  }
  double P0Tqg::Regular(double const& x) const
  {
    return 8 * CF * ( - 1 + 0.5 * x + 1 / x );
  }

  //_________________________________________________________________________________
  P0Tgq::P0Tgq():
    Expression()
  {
  }
  double P0Tgq::Regular(double const& x) const
  {
    return 2 * TR * ( 1 - 2 * x + 2 * x * x );
  }

  //_________________________________________________________________________________
  P0Tgg::P0Tgg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0Tgg::Regular(double const& x) const
  {
    return 4 * CA * ( - 2 + x - x * x + 1 / x );
  }
  double P0Tgg::Singular(double const& x) const
  {
    return 4 * CA / ( 1 - x );
  }
  double P0Tgg::Local(double const& x) const
  {
    return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
  }

  /**
   * @brief The NLO time-like splitting function classes
   */
  //_________________________________________________________________________________
  P1Tnsp::P1Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = - 40 / 9. * CF * _nf + ( 268 / 9. - 8 * zeta2 ) * CA * CF;
  }
  double P1Tnsp::Regular(double const& x) const
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
      + 4 * CF * CF * ( ( 3 * lnx / 2 + 2 * ln1mx * lnx - 2 * lnx2 ) * pqq - 5 * ( 1 - x )
			+ lnx2 * ( 1 + x ) / 2 - lnx * ( 3.5 + 3 * x / 2 ) )
      + 4 * CF * ( CF - CA / 2 ) * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    const double x1nspa = gqq1 - gqq1l;
    return x1nspa;
  }
  double P1Tnsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1Tnsp::Local(double const& x) const
  {
    const double p1delta =
      - 1 / 3. * CF * _nf + 3 / 2. * CF * CF + 17 / 6. * CA * CF + 24 * zeta3 * CF * CF - 12 * zeta3 * CA * CF
      - 8 / 3. * zeta2 * CF * _nf - 12 * zeta2 * CF * CF + 44 / 3. * zeta2 * CA * CF;
    const double x1nsc = log(1-x) * _a2 + p1delta;
    return x1nsc;
  }

  //_________________________________________________________________________________
  P1Tnsm::P1Tnsm(int const& nf):
    P1Tnsp(nf)
  {
  }
  double P1Tnsm::Regular(double const& x) const
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
      + 4 * CF * CF * ( ( 3 * lnx / 2 + 2 * ln1mx * lnx - 2 * lnx2 ) * pqq - 5 * ( 1 - x )
			+ lnx2 * ( 1 + x ) / 2 - lnx * ( 3.5 + 3 * x / 2 ) )
      - 4 * CF * ( CF - CA / 2 ) * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    const double x1nsma = gqq1 - gqq1l;
    return x1nsma;
  }

  //_________________________________________________________________________________
  P1Tps::P1Tps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1Tps::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    const double x1psa =
      _nf * CF * ( - 32 + 16 * x + 224 * x * x / 9 - 80 / 9. / x
		   - ( 20 + 36 * x ) * lnx - 32 * x * x * lnx / 3
		   + ( 4 + 4 * x ) * lnx2 );
    return x1psa;
  }

  //_________________________________________________________________________________
  P1Tqg::P1Tqg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1Tqg::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double pgq    = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
    const double pgqmx  = - ( 1 + ( 1 + x ) * ( 1 + x ) ) / x;
    const double S1x    = - dilog(1-x);
    const double S2x    = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double x1qga  =
      2 * _nf * ( 4 * CF * CF * ( - 1. / 2 + 9 * x / 2 + ( - 8 + x / 2 ) * lnx + 2 * x * ln1mx + ( 1 - x / 2 ) * lnx2
				  + ( ln1mx2 + 4 * lnx * ln1mx - 8 * S1x - 4 * Pi2 / 3 ) * pgq )
		  + 4 * CF * CA * ( 62. / 9 - 35 * x / 18 - 44 * x * x / 9 + ( 2 + 12 * x + 8 * x * x / 3 ) * lnx 
				    - 2 * x * ln1mx - ( 4 + x ) * lnx2 + pgqmx * S2x 
				    + ( - 2 * lnx * ln1mx - 3 * lnx - 3 * lnx2 / 2 - ln1mx2
					+ 8 * S1x + 7 * Pi2 / 6 + 17. / 18 ) * pgq ) );
    return x1qga;
  }

  //_________________________________________________________________________________
  P1Tgq::P1Tgq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1Tgq::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double pqg    = x * x + ( 1 - x ) * ( 1 - x );
    const double pqgmx  = x * x + ( 1 + x ) * ( 1 + x );
    const double S1x    = - dilog(1-x);
    const double S2x    = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double x1gqa  =
      + _nf * ( - 8. / 6 - ( 16. / 18 + 8 * lnx / 6 + 8 * ln1mx / 6 ) * pqg )
      + CF * ( - 2 + 3 * x + ( - 7 + 8 * x ) * lnx - 4 * ln1mx + ( 1 - 2 * x ) * lnx2 
	       + ( - 2 * pow(lnx + ln1mx, 2) - 2 * ( ln1mx - lnx ) + 16 * S1x + 2 * Pi2 - 10 ) * pqg )
      + CA * ( - 152. / 9 + 166 * x / 9 - 40 / 9. / x + ( - 4. / 3 - 76 * x / 3 ) * lnx
	       + 4 * ln1mx + ( 2 + 8 * x ) * lnx2
	       + ( 8 * lnx * ln1mx - lnx2 - 4 * lnx / 3 + 10 * ln1mx / 3
		   + 2 * ln1mx2 - 16 * S1x - 7 * Pi2 / 3 + 178. / 9 ) * pqg 
	       + 2 * pqgmx * S2x );
    return x1gqa;
  }

  //_________________________________________________________________________________
  P1Tgg::P1Tgg(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2g = - 40 / 9. * CA * _nf + ( 268 / 9. - 8 * zeta2 ) * CA * CA;
  }
  double P1Tgg::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double pgg   = ( 1 / ( 1 - x ) +  1 / x - 2 + x * ( 1 - x ) );
    const double pggmx = ( 1 / ( 1 + x ) -  1 / x - 2 - x * ( 1 + x ) );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    const double ggg1  =
      + 2 * CF * _nf * ( - 4 + 12 * x - 164 * x2 / 9 + ( 10 + 14 * x + 16 * x2 / 3 + 16 / 3. / x )
		       * lnx + 92 / 9. / x + 2 * ( 1 + x ) * lnx2 )
      + 2 * CA * _nf * ( 2 - 2 * x + 26 * ( x2 - 1 / x ) / 9 - 4 * ( 1 + x ) * lnx / 3 
			 - ( 20. / 9 + 8 * lnx / 3 ) * pgg )
      + 4 * CA * CA * ( 27 * ( 1 - x ) / 2 + 67 * ( x2 - 1 / x ) / 9
			+ ( 11. / 3 - 25 * x / 3 - 44 / 3. / x ) * lnx
			- 4 * ( 1 + x ) * lnx2
			+ ( 4 * lnx * ln1mx - 3 * lnx2 + 22 * lnx / 3
			    - 2 * zeta2 + 67. / 9 ) * pgg + 2 * pggmx * S2x );
    const double ggg1l = _a2g / ( 1 - x );
    const double x1gga = ggg1 - ggg1l;
    return x1gga;
  }
  double P1Tgg::Singular(double const& x) const
  {
    return _a2g / ( 1 - x );
  }
  double P1Tgg::Local(double const& x) const
  {
    const double p1delta = ( - 2 * CF - 8 / 3. * CA ) * _nf + ( 32 / 3. + 12 * zeta3 ) * CA * CA;
    const double x1ggc = log(1-x) * _a2g + p1delta;
    return x1ggc;
  }

  /**
   * @brief The NNLO time-like splitting function classes (parametrized)
   */
  //_________________________________________________________________________________
  P2Tnsp::P2Tnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tnsp::Regular(double const& x) const
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
      1658.7 - 707.67 * dl1 + 1327.5 * dl - 56.907 * dl * dl1 
      - 189.37 * dl_2 - 519.37 * dl1 * dl_2 - 352. / 9. * dl_3 
      + 128. / 81.* dl_4 - 4249.4 * x - 559.1 * dl1 * dl * x 
      - 1075.3 * x_2 + 593.9 * x_3
      + _nf * ( 64. / 27. * dl_3 - 176. / 81.* dl_2 - 168.89 * dl
                - 198.10 + 466.29 * x + 181.18 * x_2 - 31.84 * x_3
                + 5120. / 81. * dl1 - 50.758 * dl * dl1 + 28.551 * dl_2 * dl1
                - 39.113 * x * dl + 85.72 * x * dl * dl1 - 23.102 * x * dl_2 * dl1 )
      + _nf * _nf * ( 32.* x * dl / ( 1 - x ) * ( 3. * dl + 10. ) + 64.
		      + ( 48. * dl_2 + 352. * dl + 384. ) * ( 1 - x ) ) * d81;
    return p2nspa;
  }
  double P2Tnsp::Singular(double const& x) const
  {
    const double p2nsb = ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1 - x );
    return p2nsb;
  }
  double P2Tnsp::Local(double const& x) const
  {
    const double dl1    = log(1-x);
    const double p2nspc = 
      1174.898 * dl1 + 1295.624 + 0.001
      - _nf * ( 183.187 * dl1 + 173.938 - 0.003)
      + _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
    return p2nspc;
  }

  //_________________________________________________________________________________
  P2Tnsm::P2Tnsm(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tnsm::Regular(double const& x) const
  {
    const double x_2    = x * x;
    const double x_3    = x_2 * x;
    const double dl     = log(x);
    const double dl_2   = dl * dl;
    const double dl_3   = dl_2 * dl;
    const double dl_4   = dl_3 * dl;
    const double dl1    = log(1-x);
    const double d81    = 1. / 81.;
    const double p2nsma =
      - 140. / 81.* dl_4 - 1024. / 27. * dl_3 
      - 38.298 * dl_2 + 1625.5 * dl - 707.94 * dl1 + 1981.3
      - 4885.7 * x - 577.42 * x_2 + 407.89 * x_3
      + 1905.4 * dl_2 * dl1 + 1969.5 * x * dl_2 * dl1 
      + 4563.2 * dl * dl1 - 34.683 * x * dl_4 
      - 5140.6 * x * dl * dl1 - 437.03 * x * dl_3
      + _nf * ( 128. / 81. * dl_3 - 784. / 81. * dl_2 
		- 188.99 * dl - 217.84 + 511.92 * x + 209.19 * x_2 
		- 85.786 * x_3 + 5120. / 81. * dl1 + 71.428 * dl * dl1 
		+ 30.554 * dl_2 * dl1 + 92.453 * x * dl - 23.722 * x * dl * dl1 
		- 18.975 * x * dl_2 * dl1 )
      + _nf * _nf * ( 32. * x * dl / ( 1 - x ) * ( 3. * dl + 10. ) + 64.
		      + ( 48. * dl_2 + 352. * dl + 384. ) * ( 1 - x ) ) * d81;
    return p2nsma;
  }
  double P2Tnsm::Singular(double const& x) const
  {
    const double p2nsb = ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1 - x );
    return p2nsb;
  }
  double P2Tnsm::Local(double const& x) const
  {
    const double dl1    = log(1-x);
    const double p2nsmc =
      1174.898 * dl1 + 1295.624 - 0.002
      - _nf * ( 183.187 * dl1 + 173.938 - 0.0004 )
      + _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
    return p2nsmc;
  }

  //_________________________________________________________________________________
  P2Tnss::P2Tnss(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tnss::Regular(double const& x) const
  {
    const double x_2    = x * x;
    const double x_3    = x_2 * x;
    const double d27    = 1./27.;
    const double dl     = log(x);
    const double dl_2   = dl * dl;
    const double dl_3   = dl_2 * dl;
    const double dl_4   = dl_3 * dl;
    const double x1     = 1 - x;
    const double dl1    = log(x1);
    const double p2nssa =
      x1 * ( 151.49 + 44.51 * x - 43.12 * x_2 + 4.820 * x_3 )
      + 40. * d27 * dl_4 - 80. * d27 * dl_3 + 6.892 * dl_2
      + 178.04 * dl + dl * dl1 * ( - 173.1 + 46.18 * dl )
      + x1 * dl1 * ( - 163.9 / x - 7.208 * x );
    return _nf * p2nssa;
  }

  //_________________________________________________________________________________
  P2Tps::P2Tps(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tps::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double x_4   = x_3 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1-x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double p2ps1 =
      - 256. / ( 9. * x ) * dl_3 - 128. / ( 9. * x ) * dl_2 
      + 324.07 / x * dl + 479.87 / x
      - 5.926 * dl1_3 - 9.751 * dl1_2 - 8.65 * dl1 - 106.65
      - 848.97 * x + 368.79 * x_2 - 61.284 * x_3
      + 96.171 * dl * dl1 + 656.49 * dl + 425.14 * dl_2 
      + 47.322 * dl_3 + 9.072 * dl_4;
    const double p2ps2 =
      - 128. / ( 81. * x ) + 1.778 * dl1_2 + 16.611 * dl1 + 87.795
      - 57.688 * x - 41.827 * x_2 + 25.628 * x_3 - 7.9934 * x_4
      - 2.1031 * dl * dl1 + 57.713 * dl + 9.1682 * dl_2 
      - 1.9 * dl_3 + 0.019122 * dl_4
      + 26.294 * x * dl - 7.8645 * x * dl_3;
    const double p2psa = ( 1 - x ) * _nf * ( p2ps1 + _nf * p2ps2 );
    return p2psa;
  }

  //_________________________________________________________________________________
  P2Tqg::P2Tqg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tqg::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double x_4   = x_3 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1-x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double dl1_4 = dl1_3 * dl1;
    const double p2qg1 =
      400. / 81. * dl1_4 + 520. / 27. * dl1_3 
      - 220.13 * dl1_2 - 152.60 * dl1 + 272.85 - 7188.7 * x 
      + 5693.2 * x_2 + 146.98 * x_3 + 128.19 * x_4
      - 30.062 * dl_4 - 126.38 * dl_3 - 0.71252 * dl_2 
      + 4.4136 * dl - 1300.6 * dl * dl1 - 71.23 * dl * dl1_2 
      + 543.8 * x * dl_3 
      + 256. / x * dl_4 + 3712. / ( 3. * x ) * dl_3
      + 1001.89 / x * dl_2 + 4776.5 / x * dl + 5803.7 / x;
    const double p2qg2 =
      80. / 81.* dl1_3 + 1040./81.* dl1_2 - 16.914 * dl1
      - 871.3 + 790.13 * x - 241.23 * x_2 + 43.252 * x_3
      - 48.600 * dl_3 - 343.1 * dl_2 - 492. * dl
      + 55.048 * dl * dl1 - 4.3465 * x * dl_3
      + 6.0041 / x + 141.93 / x * dl 
      + 2912. / ( 27. * x ) * dl_2 + 1280. / ( 81. * x ) * dl_3;
    const double p2qga = 2 * _nf * ( p2qg1 + _nf * p2qg2 );
    return p2qga;
  }

  //_________________________________________________________________________________
  P2Tgq::P2Tgq(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tgq::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double x_4   = x_3 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1-x);
    const double dl1_2 = dl1 * dl1;
    const double dl1_3 = dl1_2 * dl1;
    const double dl1_4 = dl1_3 * dl1;
    const double p2gq1 =
      - 64. / x * dl_3 - 64. / x * dl_2 + 675.83/x * dl + 1141.7 / x
      + 100. / 27. * dl1_4 + 350. / 9. * dl1_3
      + 263.07 * dl1_2 + 693.84 * dl1 + 603.71
      - 882.48 * x + 4723.2 * x_2 - 4745.8 * x_3 - 175.28 * x_4
      + 1864. * dl + 1512. * dl_2 + 361.28 * dl_3 
      + 42.328 * dl_4 - 1809.4 * dl * dl1 - 107.59 * x * dl * dl1 
      - 885.5 * x * dl_4;
    const double p2gq2 =
      - 32. / ( 27. * x ) * dl_2 - 3.1752 / x * dl - 2.8986 / x
      - 100. / 27. * dl1_3 - 35.446 * dl1_2 - 103.609 * dl1
      - 113.81 + 341.26 * x - 853.35 * x_2 + 492.1 * x_3 
      + 14.803 * x_4 + 619.75 * dl + 255.62 * dl_2 
      + 21.569 * dl_3 + 966.96 * dl * dl1 - 1.593 * dl * dl1_2 
      - 333.8 * x * dl_3 - 709.1 * x * dl * dl1; 
    const double p2gq3 =
      4. / 9. * ( 4. + 6. * ( dl + dl1 )
		  + ( 1. - 2. * x + 2. * x_2 ) * ( 3.8696 + 4. * ( dl + dl1 ) 
						   + 3. * pow(dl + dl1, 2) ) );
    const double p2gqa = ( p2gq1 + _nf * p2gq2 + _nf * _nf * p2gq3 ) / 2;
    return p2gqa;
  }

  //_________________________________________________________________________________
  P2Tgg::P2Tgg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P2Tgg::Regular(double const& x) const
  {
    const double x_2   = x * x;
    const double x_3   = x_2 * x;
    const double x_4   = x_3 * x;
    const double dl    = log(x);
    const double dl_2  = dl * dl;
    const double dl_3  = dl_2 * dl;
    const double dl_4  = dl_3 * dl;
    const double dl1   = log(1-x);
    const double dl1_2 = dl1 * dl1;
    const double p2gga0 =
      + 576. / x * dl_4 + 3168. / x * dl_3 + 3651.1 / x * dl_2 
      + 10233. / x * dl + 14214.4 / x - 3590.1 * dl1 - 28489.
      + 7469. * x + 30421. * x_2 - 53017. * x_3 + 19556. * x_4
      + 191.99 * dl_4 + 3281.7 * dl_3 + 13528. * dl_2 
      + 12258. * dl - 186.4 * dl * dl1 - 21328. * dl_2 * dl1 
      + 5685.8 * x * dl_3;
    const double p2gga1 =
      + 448. / ( 9. * x ) * dl_3 + 2368. / ( 9. * x ) * dl_2 
      - 5.47 / x * dl - 804.13 / x + 248.95 + 319.97 * dl1
      + 260.6 * x + 272.79 * x_2 + 2133.2 * x_3 - 926.87 * x_4
      + 4.9934 * dl + 482.94 * dl_2 + 155.10 * dl_3 
      + 18.085 * dl_4 + 485.18 * x * dl_3 + 1266.5 * dl * dl1 
      - 29.709 * dl_2 *dl1 + 87.771 * dl * dl1_2;
    const double p2gga2 =
      + 32. / ( 27. * x ) * dl_2 + 368. / ( 81. * x ) * dl 
      + 472. / ( 243. * x )
      - 77.190 + 153.27 * x - 106.03 * x_2 + 11.995 * x_3
      - 5.0372 * dl_3 - 44.8 * dl_2 - 69.712 * dl
      - 115.01 * dl * dl1 + 96.522 * x * dl * dl1 - 62.908 * dl_2 * dl1;
    const double p2gga = p2gga0 + _nf * ( p2gga1 + _nf * p2gga2 );
    return p2gga;
  }
  double P2Tgg::Singular(double const& x) const
  {
    const double p2ggb = ( 2643.521 - _nf * 412.172 - _nf * _nf * 16. / 9. ) / ( 1 - x );
    return p2ggb;
  }
  double P2Tgg::Local(double const& x) const
  {
    const double dl1   = log(1-x);
    const double p2ggc =
      2643.521 * dl1 + 4425.448 + 0.003
      - _nf * ( 412.172 * dl1 +  528.720 - 0.001 )
      + _nf * _nf * ( - 16. / 9. * dl1 + 6.4630 - 0.0002 );
    return p2ggc;
  }

  /**
   * @brief The LO space-like splitting function for tranversely
   * polarised FFs. Reference
   * https://arxiv.org/pdf/hep-ph/0108241v1.pdf
   */
  //_________________________________________________________________________________
  P0Ttransns::P0Ttransns():
    Expression()
  {
  }
  double P0Ttransns::Regular(double const&) const
  {
    return - 4 * CF;
  }
  double P0Ttransns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double P0Ttransns::Local(double const& x) const
  {
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }

  /**
   * @brief The NLO space-like splitting function for tranversely
   * polarised FFs. Reference
   * https://arxiv.org/pdf/hep-ph/0108241v1.pdf
   */
  //_________________________________________________________________________________
  P1Ttransnsp::P1Ttransnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = 2 * CA * CF * ( 134. / 9. - 2 * Pi2 / 3 ) - 80 * _nf * CF * TR / 9;
  }
  double P1Ttransnsp::Regular(double const& x) const
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
  double P1Ttransnsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1Ttransnsp::Local(double const& x) const
  {
    const double p1delta =
      + 4 * CF * CF * ( 3. / 8. - Pi2 / 2 + 6 * zeta3 )
      + 2 * CF * CA * ( 17. / 12. + 11 * Pi2 / 9 - 6 * zeta3 )
      - 8 * _nf * CF * TR * ( 1. / 4. + Pi2 / 3 );
    const double x1nsc = log(1-x) * _a2 + p1delta;
    return x1nsc;
  }

  //_________________________________________________________________________________
  P1Ttransnsm::P1Ttransnsm(int const& nf):
    P1Ttransnsp(nf)
  {
  }
  double P1Ttransnsm::Regular(double const& x) const
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
