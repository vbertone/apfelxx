//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"
#include "apfel/tools.h"
#include "apfel/specialfunctions.h"

#include <cmath>

using namespace std;

namespace apfel
{
  /**
   * @brief The LO splitting function classes
   */
  class P0ns: public Expression
  {
  public:
  P0ns(): Expression() { }
    double Regular(double const& x)  const { return - 2 * CF * ( 1 + x ); }
    double Singular(double const& x) const { return 4 * CF / ( 1 - x ); }
    double Local(double const& x)    const { return 4 * CF * log( 1 - x ) + 3 * CF; }
  };

  class P0qg: public Expression
  {
  public:
  P0qg(): Expression() { }
    double Regular(double const& x)  const { return 2 * ( 1 - 2 * x + 2 * x * x ); }
  };

  class P0gq: public Expression
  {
  public:
  P0gq(): Expression() { }
    double Regular(double const& x)  const { return 4 * CF * ( - 1 + 0.5 * x + 1 / x ); }
  };

  class P0gg: public Expression
  {
  public:
  P0gg(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x)  const { return 4 * CA * ( - 2 + x - x * x + 1 / x ); }
    double Singular(double const& x) const { return 4 * CA / ( 1 - x ); }
    double Local(double const& x)    const { return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA; }
  private:
    int const _nf;
  };

  /**
   * @brief The NLO splitting function classes
   */
  class P1nsp: public Expression
  {
  public:
  P1nsp(int const& nf): Expression(), _nf(nf) { _a2 = - 40 / 9. * CF * _nf + ( 268 / 9. - 8 * zeta2 ) * CA * CF; }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pqq   = 2 / ( 1 - x ) - 1 - x;
      const auto pqqmx = 2 / ( 1 + x ) - 1 + x;
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
      const auto gqq1  =
	+ 2 * CF * _nf * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 ) 
	+ 4 * CA * CF * ( ( 67 / 18. + 11 * lnx / 6 + lnx2 / 2 - Pi2 / 6 ) * pqq
			  + 20 * ( 1 - x ) / 3 + lnx * ( 1 + x ) )
	+ 4 * CF * CF * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x )
			  - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
	+ 4 * CF * ( CF - CA / 2 ) * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
      const auto gqq1l = _a2 / ( 1 - x );
      const auto x1nspa = gqq1 - gqq1l;
      return x1nspa;
    }
    virtual double Singular(double const& x) const
    {
      return _a2 / ( 1 - x );
    }
    virtual double Local(double const& x) const
    {
      const auto p1delta =
	- 1 / 3. * CF * _nf + 3 / 2. * CF * CF + 17 / 6. * CA * CF + 24 * zeta3 * CF * CF - 12 * zeta3 * CA * CF
	- 8 / 3. * zeta2 * CF * _nf - 12 * zeta2 * CF * CF + 44 / 3. * zeta2 * CA * CF;
      const auto x1nsc = log(1-x) * _a2 + p1delta;
      return x1nsc;
    }
  protected:
    int const _nf;
    double    _a2;
  };

  class P1nsm: public P1nsp
  {
  public:
  P1nsm(int const& nf): P1nsp(nf) { }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pqq   = 2 / ( 1 - x ) - 1 - x;
      const auto pqqmx = 2 / ( 1 + x ) - 1 + x;
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
      const auto gqq1  =
	+ 2 * CF * _nf * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 ) 
	+ 4 * CA * CF * ( ( 67 / 18. + 11 * lnx / 6 + lnx2 / 2 - Pi2 / 6 ) * pqq
			  + 20 * ( 1 - x ) / 3 + lnx * ( 1 + x ) )
	+ 4 * CF * CF * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x )
			  - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
	- 4 * CF * ( CF - CA / 2 ) * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
      const auto gqq1l = _a2 / ( 1 - x );
      const auto x1nsma = gqq1 - gqq1l;
      return x1nsma;
    }
  };

  class P1ps: public Expression
  {
  public:
  P1ps(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto lnx  = log(x);
      const auto lnx2 = lnx * lnx;
      const auto x1psa =
	_nf * CF * ( - 8 + 24 * x - 224 / 9. * x * x + 80 / 9. / x
		     + ( 4 + 20 * x ) * lnx + 32 / 3. * x * x * lnx
		     - ( 4 + 4 * x ) * lnx2 );
      return x1psa;
    }
  private:
    int const _nf;
  };

  class P1qg: public Expression
  {
  public:
  P1qg(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pqg   = x * x + ( 1 - x ) * ( 1 - x );
      const auto pqgmx = x * x + ( 1 + x ) * ( 1 + x );
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
      const auto x1qga =
	+ 2 * CF * _nf * ( 4  + 4 * ln1mx + ( 10 - 4 * ( ln1mx - lnx ) + 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 2 * Pi2 / 3 ) * pqg
			   - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) - 9 * x )
	+ 2 * CA * _nf * ( 182 / 9. - 4 * ln1mx
			   + ( - 218 / 9. + 4 * ln1mx - 2 * ln1mx * ln1mx + 44 * lnx / 3 - lnx2 + Pi2 / 3 ) * pqg
			   + 2 * pqgmx * S2x + 40 / ( 9 * x ) + 14 * x / 9 - lnx2 * ( 2 + 8 * x )
			   + lnx * ( - 38 / 3. + 136 * x / 3 ) );
      return x1qga;
    }
  private:
    int const _nf;
  };

  class P1gq: public Expression
  {
  public:
  P1gq(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pgq   = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
      const auto pgqmx = - ( 1 + ( 1 + x ) * ( 1 + x ) ) / x;
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
      const auto x1gqa =
	+ 2 * CF * _nf * ( - ( 20 / 9. + 4 * ln1mx / 3 ) * pgq - 4 * x / 3 )
	+ 4 * CF * CF  * ( - 2.5 - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq - lnx2 * ( 1 - x / 2 ) - 7 * x / 2 
			   - 2 * ln1mx * x + lnx * ( 2 + 7 * x / 2 ) )
	+ 4 * CA * CF  * ( 28 / 9. + pgq * ( 0.5 + 11 * ln1mx / 3 + ln1mx * ln1mx - 2 * ln1mx * lnx + lnx2 / 2 - Pi2 / 6 ) + pgqmx * S2x 
			   + 65 * x / 18 + 2 * ln1mx * x + 44 * x * x / 9 + lnx2 * ( 4 + x ) - lnx * ( 12 + 5 * x + 8 * x * x / 3 ) );
      return x1gqa;
    }
  private:
    int const _nf;
  };

  class P1gg: public Expression
  {
  public:
  P1gg(int const& nf): Expression(), _nf(nf) { _a2g = - 40 / 9. * CA * _nf + ( 268 / 9. - 8 * zeta2 ) * CA * CA; }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pgg   = ( 1 / ( 1 - x ) +  1 / x - 2 + x * ( 1 - x ) );
      const auto pggmx = ( 1 / ( 1 + x ) -  1 / x - 2 - x * ( 1 + x ) );
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
      const auto ggg1  =
	+ 2 * CF * _nf * ( - 16 + 4 / ( 3 * x ) + 8 * x + ( 20 * x * x ) / 3 - lnx2 * ( 2 + 2 * x ) - lnx * ( 6 + 10 * x ) )
	+ 2 * CA * _nf * ( 2 - 20 * pgg / 9 - 2 * x - 4 * lnx * ( 1 + x ) / 3 + 26 * ( - 1 / x + x * x ) / 9 )
	+ 4 * CA *  CA * ( pgg * ( 67 / 9. - 4 * ln1mx * lnx + lnx2 - Pi2 / 3 ) + 2 * pggmx * S2x
			+ 27 * ( 1 - x ) / 2 + 4 * lnx2 * ( 1 + x ) + 67 * ( - 1 / x + x * x ) / 9
			- lnx * ( 25 / 3. - 11 * x / 3 + 44 * x * x / 3 ) );
      const auto ggg1l = _a2g / ( 1 - x );
      const auto x1gga = ggg1 - ggg1l;
      return x1gga;
    }
    double Singular(double const& x) const
    {
      return _a2g / ( 1 - x );
    }
    double Local(double const& x) const
    {
      const auto p1delta = ( - 2 * CF - 8 / 3. * CA ) * _nf + ( 32 / 3. + 12 * zeta3 ) * CA * CA;
      const auto x1ggc = log(1-x) * _a2g + p1delta;
      return x1ggc;
    }
  private:
    int const _nf;
    double    _a2g;
  };

  /**
   * @brief The NNLO splitting function classes (parametrized)
   */
  class P2nsp: public Expression
  {
  public:
  P2nsp(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto x_2    = pow(x, 2);
      const auto x_3    = x_2 * x;
      const auto dl     = log(x);
      const auto dl_2   = pow(dl,2);
      const auto dl_3   = dl_2 * dl;
      const auto dl_4   = dl_3 * dl;
      const auto dl1    = log(1-x);
      const auto d81    = 1. / 81.;
      const auto p2nspa =
        1641.1 - 3135. * x + 243.6 * x_2 - 522.1 * x_3
        + 128. * d81 * dl_4 + 2400. * d81 * dl_3
        + 294.9 * dl_2 + 1258. * dl
	+ 714.1 * dl1 + dl * dl1 * ( 563.9 + 256.8 * dl )
	+ _nf * ( -197.0 + 381.1 * x + 72.94 * x_2 + 44.79 * x_3
		 - 192. * d81 * dl_3  - 2608. * d81 * dl_2 - 152.6 * dl
		 - 5120. * d81 * dl1 - 56.66 * dl * dl1 - 1.497 * x * dl_3 )
	+ _nf * _nf * ( 32. * x * dl / ( 1. - x ) * ( 3. * dl + 10. ) + 64.
			+ ( 48. * dl_2 + 352. * dl + 384. ) * ( 1. - x ) ) * d81;
      return p2nspa;
    }
    double Singular(double const& x) const
    {
      const auto p2nsb = ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1. - x );
      return p2nsb;
    }
    double Local(double const& x) const
    {
      const auto dl1    = log(1.-x);
      const auto p2nspc =
	1174.898 * dl1 + 1295.624 - 0.24
	- _nf * ( 183.187 * dl1 + 173.938 - 0.011 )
	+ _nf * _nf * ( - 64. / 81. * dl1 + 1.13067 );
      return p2nspc;
    }
  private:
    int const _nf;
  };

  class P2nsm: public Expression
  {
  public:
  P2nsm(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
       const auto x_2    = pow(x, 2);
       const auto x_3    = x_2 * x;
       const auto dl     = log(x);
       const auto dl_2   = pow(dl,2);
       const auto dl_3   = dl_2 * dl;
       const auto dl1    = log(1.-x);
       const auto d81    = 1./81.;       
       const auto p2nsma =
         1860.2 - 3505.* x + 297.0 * x_2 - 433.2 * x_3
         + 116. * d81 * dl_3 * dl + 2880. * d81 * dl_3
         + 399.2 * dl_2 + 1465.2 * dl
	 + 714.1 * dl1 + dl * dl1 * ( 684.0 + 251.2 * dl )
	 + _nf * ( -216.62 + 406.5 * x + 77.89 * x_2 + 34.76 * x_3
		   - 256. * d81 * dl_3  - 3216. * d81 * dl_2 - 172.69 * dl
		   - 5120. * d81 * dl1 - 65.43 * dl * dl1 - 1.136 * x * dl_3 )
	 + _nf * _nf * ( 32.* x * dl / (1.-x) * ( 3. * dl + 10. ) + 64.
			 + ( 48.* dl_2 + 352.* dl + 384. ) * ( 1.-x ) ) * d81;
      return p2nsma;
    }
    double Singular(double const& x) const
    {
      const auto p2nsb = ( 1174.898 - _nf * 183.187 - _nf * _nf * 64. / 81. ) / ( 1. - x );
      return p2nsb;
    }
    double Local(double const& x) const
    {
       const auto dl1    = log(1.-x);
       const auto p2nsmc =
	 1174.898 * dl1 + 1295.624 - 0.154
	 - _nf * ( 183.187 * dl1 + 173.938  - 0.005 )
	 + _nf * _nf * ( - 64./81. * dl1 + 1.13067 );
      return p2nsmc;
    }
  private:
    int const _nf;
  };

  class P2nss: public Expression
  {
  public:
  P2nss(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto x_2    = pow(x,2);
      const auto d27    = 1./27.;
      const auto dl     = log(x);
      const auto dl_2   = pow(dl,2);
      const auto dl_3   = dl_2 * dl;
      const auto dl_4   = dl_3 * dl;
      const auto x1     = 1.- x;
      const auto dl1    = log(x1);
      const auto p2nssa =
        x1 * ( 151.49 + 44.51 * x - 43.12 * x_2 + 4.820 * x_2 * x )
        + 40. * d27 * dl_4 - 80. * d27 * dl_3 + 6.892 * dl_2
	+ 178.04 * dl + dl * dl1 * ( - 173.1 + 46.18 * dl )
	+ x1 * dl1 * ( - 163.9 / x - 7.208 * x );
      return _nf * p2nssa;
    }
  private:
    int const _nf;
  };

  class P2ps: public Expression
  {
  public:
  P2ps(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto x_2   = pow(x,2);
      const auto x_3   = x_2 * x;

      const auto dl    = log(x);
      const auto dl_2  = pow(dl, 2);
      const auto dl_3  = dl_2 * dl;
      const auto dl_4  = dl_3 * dl;

      const auto dl1   = log(1.-x);
      const auto dl1_2 = pow(dl1, 2);
      const auto dl1_3 = dl1_2 * dl1;

      const auto  p2ps1 =
        - 3584. / ( 27. * x ) * dl - 506.0 / x + 160. / 27. * dl_4
        - 400. / 9. * dl_3 + 131.4 * dl_2 - 661.6 * dl
        - 5.926  * dl1_3 - 9.751 * dl1_2 - 72.11 * dl1
        + 177.4 + 392.9 * x - 101.4 * x_2 - 57.04 * dl * dl1;
      const auto p2ps2  =
        256. / ( 81. * x ) + 32. / 27. * dl_3 + 17.89 * dl_2
        + 61.75 * dl + 1.778 * dl1_2 + 5.944 * dl1 + 100.1
        - 125.2 * x + 49.26 * x_2 - 12.59 * x_3
	- 1.889 * dl * dl1;
      const auto p2psa = ( 1. - x ) * _nf * ( p2ps1 + _nf * p2ps2 );
      return p2psa;
    }
  private:
    int const _nf;
  };

  class P2qg: public Expression
  {
  public:
  P2qg(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto x_2   = pow(x,2);
      const auto x_3   = x_2 * x;

      const auto dl    = log(x);
      const auto dl_2  = pow(dl, 2);
      const auto dl_3  = dl_2 * dl;
      const auto dl_4  = dl_3 * dl;

      const auto dl1   = log(1.-x);
      const auto dl1_2 = pow(dl1, 2);
      const auto dl1_3 = dl1_2 * dl1;
      const auto dl1_4 = dl1_3 * dl1;

      const auto p2qg1 =
        - 896. / ( 3. * x ) * dl - 1268.3 / x + 536./27. * dl_4
        - 44. / 3. * dl_3 + 881.5 * dl_2 + 424.9 * dl
        + 100. / 27. * dl1_4 - 70. / 9. * dl1_3
        - 120.5 * dl1_2 + 104.42 * dl1
        + 2522. - 3316. * x + 2126. * x_2
        + dl * dl1 * ( 1823. - 25.22 * dl ) - 252.5 * x * dl_3;
      const auto p2qg2 =
        1112. / ( 243. * x ) - 16. / 9. * dl_4
        - 376. / 27. * dl_3 - 90.8 * dl_2 - 254.0 * dl
        + 20./27. * dl1_3 + 200. / 27. * dl1_2 - 5.496 * dl1
        - 252.0  + 158.0 * x + 145.4 * x_2 - 139.28 * x_3
        - dl * dl1 * ( 53.09  + 80.616 * dl ) - 98.07 * x * dl_2
        + 11.70 * x * dl_3;
      const auto p2qga = _nf * ( p2qg1 + _nf * p2qg2 );
      return p2qga;      
    }
  private:
    int const _nf;
  };

  class P2gq: public Expression
  {
  public:
  P2gq(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto x_2   = pow(x,2);
      const auto x_3   = x_2 * x;

      const auto dl    = log(x);
      const auto dl_2  = pow(dl, 2);
      const auto dl_3  = dl_2 * dl;
      const auto dl_4  = dl_3 * dl;

      const auto dl1   = log(1.-x);
      const auto dl1_2 = pow(dl1, 2);
      const auto dl1_3 = dl1_2 * dl1;
      const auto dl1_4 = dl1_3 * dl1;

      const auto p2gq0 =
        1189.3 * dl / x + 6163.1 / x - 4288. / 81. * dl_4
        + 1568. / 9. * dl_3 - 1794. * dl_2 + 4033. * dl
        + 400. / 81. * dl1_4 + 2200. / 27. * dl1_3
        + 606.3 * dl1_2 + 2193. * dl1
        - 4307. + 489.3 * x + 1452.* x_2 + 146.0 * x_3
        - 447.3 * dl_2 * dl1 - 972.9 * x * dl_2;
     const auto p2gq1 =
       71.082 * dl / x  - 46.41 / x + 128. / 27. * dl_4
       + 704/81. * dl_3 + 20.39 * dl_2 + 174.8 * dl
       - 400./81. * dl1_3 - 68.069 * dl1_2 - 296.7 * dl1
       - 183.8 + 33.35 * x - 277.9 * x * x + 108.6 * x * dl_2
       - 49.68 * dl * dl1;
     const auto p2gq2 =
       ( 64. * ( - 1. / x + 1. + 2. * x )
	 + 320. * dl1 * ( 1. / x - 1. + 0.8 * x )
	 + 96. * dl1_2 * ( 1. / x - 1. + 0.5 * x ) ) / 27.;
     const auto p2gqa = ( p2gq0 + _nf * ( p2gq1 + _nf * p2gq2 ) );
     return p2gqa;
    }
  private:
    int const _nf;
  };

  class P2gg: public Expression
  {
  public:
  P2gg(int const& nf): Expression(), _nf(nf) { }
    double Regular(double const& x) const
    {
      const auto x_2   = pow(x,2);
      const auto x_3   = x_2 * x;

      const auto dl    = log(x);
      const auto dl_2  = pow(dl, 2);
      const auto dl_3  = dl_2 * dl;
      const auto dl_4  = dl_3 * dl;

      const auto dl1   = log(1.-x);

      const auto p2gga0 =
        2675.8 * dl / x + 14214. / x - 144. * dl_4 + 72. * dl_3
        - 7471. * dl_2 + 274.4 * dl + 3589. * dl1 - 20852.
        + 3968. * x - 3363. * x_2 + 4848. * x_3
	+ dl * dl1 * ( 7305. + 8757. * dl );
      const auto p2gga1 =
        157.27 * dl / x + 182.96 / x + 512./27. * dl_4
        + 832. / 9. * dl_3 + 491.3 * dl_2 + 1541. * dl
        - 320.0 * dl1 - 350.2 + 755.7 * x - 713.8 * x_2
        + 559.3 * x_3 + dl * dl1 * ( 26.15 - 808.7 * dl );
     const auto p2gga2 =
       - 680. / ( 243. * x ) - 32. / 27. * dl_3 + 9.680 * dl_2
       - 3.422 * dl - 13.878 + 153.4 * x - 187.7 * x_2
       + 52.75 * x_3 - dl * dl1 * ( 115.6 - 85.25 * x + 63.23 * dl);
      const auto p2gga = p2gga0 + _nf * ( p2gga1 + _nf * p2gga2 );
      return p2gga;
    }
    double Singular(double const& x) const
    {
      const auto p2ggb = ( 2643.521 - _nf * 412.172 - _nf * _nf * 16. / 9. ) / ( 1.- x );
      return p2ggb;
    }
    double Local(double const& x) const
    {
      const auto dl1   = log(1.-x);
      const auto p2ggc =
	2643.521 * dl1 + 4425.448 + 0.446
	- _nf * ( 412.172 * dl1 +  528.720 + 0.003 )
	+ _nf * _nf * ( - 16. / 9. * dl1 + 6.4630 );
      return p2ggc;
    }
  private:
    int const _nf;
  };

}
