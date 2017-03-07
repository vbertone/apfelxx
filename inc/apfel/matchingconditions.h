//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"

using namespace std;

namespace apfel
{
  /**
   * @brief Identity expression (delta function)
   */
  class Identity: public Expression
  {
  public:
  Identity(): Expression() { }
    double Local(double const& x) const { return 1 + 0 * x; }
  };

  /**
   * @brief Zero expression
   */
  class Null: public Expression
  {
  public:
  Null(): Expression() { }
  };

  /**
   * @notes Expressions taken from https://arxiv.org/pdf/hep-ph/9612398.pdf.
   * @notes To note that in these expressions ln(m2/mu2) appears while we need ln(mu2/m2),
   * @notes so we need to include a sign minus in front of every term linear in this log.
   */

  /**
   * @brief O(alpha_s) mathcing conditions
   */
  // Term propotional to ln(mu2/m2) of eq. (B.2)
  class AS1Hg_L: public Expression
  {
  public:
  AS1Hg_L(): Expression() { }
    double Regular(double const& x) const { return  4 * TR * ( x * x + ( 1 - x ) * ( 1 - x ) ); }
  };

  // Term propotional to ln(mu2/m2) of eq (B.6)
  class AS1ggH_L: public Expression
  {
  public:
  AS1ggH_L(): Expression() { }
    double Local(double const& x) const { return  - 4 * TR / 3. + 0 * x; }
  };

  /**
   * @brief O(alpha_s^2) mathcing conditions
   */
  // Constant term of eq (B.1)
  class APS2Hq_0: public Expression
  {
  public:
  APS2Hq_0(): Expression() { }
    double Regular(double const& x) const
    {
      const auto lnx    = log(x);
      const auto x2     = x * x;
      const auto Li21mx = dilog(1-x);
      const auto S121mx = wgplg(1,2,1-x);
      const auto a0 =
	( 1 + x ) * ( 32 * S121mx + 16 * lnx * Li21mx - 16 * zeta2 * lnx - 4 * pow(lnx,3) / 3 ) 
        + ( 32 / 3. / x + 8 - 8 * x - 32 * x2 / 3 ) * ( Li21mx - zeta2 )
	+ ( 2 + 10 * x + 16 * x2 / 3 ) * lnx * lnx - ( 56 / 3. + 88 * x / 3 + 448 * x2 / 9 ) * lnx
	- 448 / 27. / x - 4 / 3. - 124 * x / 3 + 1600 * x2 / 27;
      return CF * TR * a0;
    }
  };

  // Constant term of eq (B.3)
  class AS2Hg_0: public Expression
  {
  public:
  AS2Hg_0(): Expression() { }
    double Regular(double const& x) const
    {
      const auto S121mx = wgplg(1,2,1-x);
      const auto S12mx  = wgplg(1,2,-x);
      const auto S211mx = wgplg(2,1,1-x);
      const auto S21mx  = wgplg(2,1,-x);
      const auto S111mx = dilog(1-x);
      const auto S11mx  = dilog(-x);

      const auto x2     = x * x;
      const auto lnx    = log(x);
      const auto lnx2   = lnx * lnx;
      const auto lnx3   = lnx2 * lnx;
      const auto ln1mx  = log(1-x);
      const auto ln1mx2 = ln1mx * ln1mx; 
      const auto ln1mx3 = ln1mx2 * ln1mx; 
      const auto ln1px  = log(1+x);
      const auto ln1px2 = ln1px * ln1px;

      // CF * TR  part
      const auto a01 =
	( 1 - 2 * x + 2 * x2 ) * ( 8 * zeta3 + 4 * ln1mx3 / 3 - 8 * ln1mx * S111mx
				   + 8 * zeta2 * lnx - 4 * lnx * ln1mx2 + 2 * lnx3 / 3 
				   - 8 * lnx * S111mx + 8 * S211mx - 24 * S121mx );
      const auto b01 =
	- ( 4 + 96 * x - 64 * x2 ) * S111mx 
	- ( 4 - 48 * x + 40 * x2 ) * zeta2 
	- ( 8 + 48 * x - 24 * x2 ) * lnx * ln1mx 
	+ ( 4 + 8 * x - 12 * x2 ) * ln1mx2 
	- ( 1 + 12 * x - 20 * x2 ) * lnx2 
	- ( 52 * x - 48 * x2 ) * ln1mx 
	- ( 16 + 18 * x + 48 * x2 ) * lnx 
	+ 26 - 82 * x + 80 * x2
	+ x2 * ( - 16 * zeta2 * lnx + 4 * lnx3 / 3 +  16 * lnx * S111mx +  32 * S121mx );
      // CA * TR  part
      const auto a02 =
	( 1 - 2 * x + 2 * x2 ) * ( - 4 * ln1mx3 / 3 + 8 * ln1mx * S111mx - 8 * S211mx ) 
	+ ( 1 + 2 * x + 2 * x2 ) * ( - 8 * zeta2 * ln1px - 16 * ln1px * S11mx - 8 * lnx * ln1px2 
				     + 4 * lnx2 * ln1px + 8 * lnx * S11mx - 8 * S21mx - 16 * S12mx )
	+ ( 16 + 64 * x ) * ( 2 * S121mx + lnx * S111mx )
	- ( 4 + 8 * x ) * lnx3 / 3 
	+ ( 8 - 32 * x + 16 * x2 ) * zeta3 
	- ( 16 + 64 * x ) * zeta2 * lnx;
      const auto b02 =
	( 16 * x + 16 * x2 ) * ( S11mx + lnx * ln1px ) // There is a typo here in the e-Print ( 16 -> 16x )
	+ ( 32 / x / 3 + 12 + 64 * x - 272 * x2 / 3 ) * S111mx
	- ( 12 + 48 * x - 260 * x2 / 3 + 32 / x / 3 ) * zeta2
	- 4 * x2 * lnx * ln1mx 
	- ( 2 + 8 * x - 10 * x2 ) * ln1mx2 
	+ ( 2 + 8 * x + 46 * x2 / 3 ) * lnx2 
	+ ( 4 + 16 * x - 16 * x2 ) * ln1mx
	- ( 56. / 3 + 172 * x / 3 + 1600 * x2 / 9 ) * lnx 
	- 448 / x / 27 - 4. / 3 - 628 * x / 3 + 6352 * x2 / 27;

      return TR * ( CF * ( a01 + b01 )  +  CA * ( a02 + b02 ) );
    }
  };

  // Constant term of eq (B.4)
  class ANS2qqH_0: public Expression
  {
  public:
  ANS2qqH_0(): Expression() { }
    double Regular(double const& x) const
    {
      const auto x2    = x * x;
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto a0    =
	( 1 + x2 ) * ( 2 * lnx2 / 3 + 20 * lnx / 9 ) / ( 1 - x ) 
        + 8 * ( 1 - x ) * lnx / 3 + 44 / 27. - 268 * x / 27;
      return CF * TR * a0;
    }
    double Singular(double const& x) const
    {
      const auto a0 = 224 / 27.;
      return CF * TR * a0 / ( 1 - x );
    }
    double Local(double const& x) const
    {
      const auto ln1mx = log(1-x);
      const auto a0 = - 8 * zeta3 / 3 + 40 * zeta2 / 9 + 73 / 18. + 224 * ln1mx / 27;
      return CF * TR * a0;
    }
  };

  // Constant term of eq (B.5)
  class AS2gqH_0: public Expression
  {
  public:
  AS2gqH_0(): Expression() { }
    double Regular(double const& x) const
    {
      const auto ln1mx = log(1 - x);
      const auto a0    =
	4 * ( 2 / x - 2 + x ) * ln1mx * ln1mx / 3 
        + 8 * ( 10 / x - 10 + 8 * x ) * ln1mx / 9 
        + ( 448 / x - 448 + 344 * x ) / 27;
      return CF * TR * a0;
    }
  };

  // Constant term of eq (B.7)
  class AS2ggH_0: public Expression
  {
  public:
  AS2ggH_0(): Expression() { }
    double Regular(double const& x) const
    {
      const auto x2    = x * x;
      const auto lnx   = log(x);
      const auto lnx2  = lnx * lnx;
      const auto lnx3  = lnx2 * lnx;
      const auto ln1mx = log(1-x);
      // CF * TR part
      const auto a01   =
	4 * ( 1 + x ) * lnx3 / 3 + ( 6 + 10 * x) * lnx2 
	+ ( 32 + 48 * x ) * lnx - 8 / x + 80 - 48 * x - 24 * x2;
      // CA * TR part
      const auto a02   =
	4 * ( 1 + x ) * lnx2 / 3 + ( 52 + 88 * x ) * lnx / 9
	- 4 * x * ln1mx / 3 + ( 556 / x - 628 + 548 * x - 700 * x2 ) / 27; 
      return TR * ( CF * a01 + CA * a02 );
    }
    double Singular(double const& x) const
    {
      const auto a0 = 224 / 27.;
      return CA * TR * a0 / ( 1 - x );
    }
    double Local(double const& x) const
    {
      const auto ln1mx = log(1-x);
      const auto a01 = - 15.;
      const auto a02 = 10 / 9. + 224 * ln1mx / 27;
      return TR * ( CF * a01 + CA * a02 );
    }
  };


}