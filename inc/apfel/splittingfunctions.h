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
  P1nsp(int const& nf): Expression(), _nf(nf) { _a2 = - 40 / 9 * CF * _nf + ( 268 / 9 - 8 * zeta2 ) * CA * CF; }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pqq   = 2 / ( 1 - x ) - 1 - x;
      const auto pqqmx = 2 / ( 1 + x ) - 1 + x;
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - M_PI / 6;
      const auto gqq1 =
	+ 2 * CF * _nf * ( ( - 1.1111111111111112 - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 ) 
	+ 4 * CA * CF * ( ( 3.7222222222222223 + 11 * lnx / 6 + lnx2 / 2 - Pi2 / 6 ) * pqq
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
	- 1 / 3 * CF * _nf + 3 / 2 * CF * CF + 17 /6 * CA * CF + 24 * zeta3 * CF * CF - 12 * zeta3 * CA * CF
	- 8 /3 * zeta2 * CF * _nf - 12 * zeta2 * CF * CF + 44 / 3 * zeta2 * CA * CF;
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
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - M_PI / 6;
      const auto gqq1 = 
	+ 2 * CF * _nf * ( ( - 1.1111111111111112 - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 ) 
	+ 4 * CA * CF * ( ( 3.7222222222222223 + 11 * lnx / 6 + lnx2 / 2 - Pi2 / 6 ) * pqq
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
	_nf * CF * ( - 8 + 24 * x - 224 / 9 * x * x + 80 / 9 / x
		     + ( 4 + 20 * x ) * lnx + 32 / 3 * x * x * lnx
		     - ( 4 - 4 * x ) * lnx2 );
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
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - M_PI / 6;
      const auto x1qga =
	+ 2 * CF * _nf * ( 4  + 4 * ln1mx + ( 10 - 4 * ( ln1mx - lnx ) + 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 2 * Pi2 / 3 ) * pqg
			   - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) - 9 * x )
	+ 2 * CA * _nf * ( 20.22222222222222 - 4 * ln1mx
			   + ( - 24.22222222222222 + 4 * ln1mx - 2 * ln1mx * ln1mx + 44 * lnx /3 - lnx2 + Pi2 / 3 ) * pqg
			   + 2 * pqgmx * S2x + 40 / ( 9 * x ) + 14 * x / 9 - lnx2 * ( 2 + 8 * x )
			   + lnx * ( - 12.666666666666666 + 136 * x / 3 ) );
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
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - M_PI / 6;
      const auto x1gqa = 
	+ 2 * CF * _nf * ( - ( 2.2222222222222223 + 4 * ln1mx / 3 ) * pgq - 4 * x / 3 )
	+ 4 * CF * CF  * ( - 2.5 - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq - lnx2 * ( 1 - x / 2 ) - 7 * x / 2 
			   - 2 * ln1mx * x + lnx * ( 2 + 7 * x / 2 ) )
	+ 4 * CA * CF  * ( 3.111111111111111
			   + pgq * ( 0.5 + 11 * ln1mx / 3 + ln1mx * ln1mx - 2 * ln1mx * lnx + lnx2 / 2 - Pi2 / 6 ) + pgqmx * S2x 
			   + 65 * x / 18 + 2 * ln1mx * x + 44 * x * x / 9 + lnx2 * ( 4 + x ) - lnx * ( 12 + 5 * x + 8 * x * x / 3 ) );
      return x1gqa;
    }
  private:
    int const _nf;
  };

  class P1gg: public Expression
  {
  public:
  P1gg(int const& nf): Expression(), _nf(nf) { _a2g = - 40 / 9 * CA * _nf + ( 268 / 9 - 8 * zeta2 ) * CA * CA; }
    double Regular(double const& x) const
    {
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto ln1mx = log(1-x);
      const auto pgg   = ( 1 / ( 1 - x ) +  1 / x - 2 + x * ( 1 - x ) );
      const auto pggmx = ( 1 / ( 1 + x ) -  1 / x - 2 - x * ( 1 + x ) );
      const auto S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - M_PI / 6;
      const auto ggg1  =
	+ 2 * CF * _nf * ( - 16 + 4 / ( 3 * x ) + 8 * x + ( 20 * x * x ) / 3 - lnx * x * ( 2 + 2 * x ) - lnx * ( 6 + 10 * x ) )
	+ 2 * CA * _nf * ( 2 - ( 20 * pgg ) / 9 - 2 * x - ( 4 * lnx * ( 1 + x ) ) / 3 + ( 26 * ( - ( 1 / x ) + x * x ) ) / 9 )
	+ 4 * CA * CA * ( pgg * ( 7.444444444444445 - 4 * ln1mx * lnx + lnx * x - Pi2 / 3 ) + 2 * pggmx * S2x
			  + ( 27 * ( 1 - x ) ) / 2 + 4 * lnx * x * ( 1 + x ) + ( 67 * ( - ( 1 / x ) + x * x ) ) / 9
			  - lnx * ( 8.333333333333334 - ( 11 * x ) / 3 + ( 44 * x * x ) / 3 ) );
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
      const auto p1delta = ( - 2 * CF - 8 / 3 * CA ) * _nf + ( 32 / 3 + 12 * zeta3 ) * CA * CA;
      const auto x1ggc = log(1-x) * _a2g + p1delta;
      return x1ggc;
    }
  private:
    int const _nf;
    double    _a2g;
  };

}
