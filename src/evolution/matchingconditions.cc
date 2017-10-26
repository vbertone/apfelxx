//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/matchingconditions.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  AS1Hg_L::AS1Hg_L():
    Expression()
  {
  }
  double AS1Hg_L::Regular(double const& x) const
  {
    return 4 * TR * ( x * x + ( 1 - x ) * ( 1 - x ) );
  }

  //_________________________________________________________________________________
  AS1ggH_L::AS1ggH_L():
    Expression()
  {
  }
  double AS1ggH_L::Local(double const&) const
  {
    return - 4 * TR / 3.;
  }

  //_________________________________________________________________________________
  APS2Hq_0::APS2Hq_0():
    Expression()
  {
  }
  double APS2Hq_0::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double x2     = x * x;
    const double Li21mx = dilog(1-x);
    const double S121mx = wgplg(1,2,1-x);
    const double a0 =
      ( 1 + x ) * ( 32 * S121mx + 16 * lnx * Li21mx - 16 * zeta2 * lnx - 4 * pow(lnx,3) / 3 ) 
      + ( 32 / 3. / x + 8 - 8 * x - 32 * x2 / 3 ) * ( Li21mx - zeta2 )
      + ( 2 + 10 * x + 16 * x2 / 3 ) * lnx * lnx - ( 56 / 3. + 88 * x / 3 + 448 * x2 / 9 ) * lnx
      - 448 / 27. / x - 4 / 3. - 124 * x / 3 + 1600 * x2 / 27;
    return CF * TR * a0;
  }

  //_________________________________________________________________________________
  AS2Hg_0::AS2Hg_0():
    Expression()
  {
  }
  double AS2Hg_0::Regular(double const& x) const
  {
    const double S121mx = wgplg(1,2,1-x);
    const double S12mx  = wgplg(1,2,-x);
    const double S211mx = wgplg(2,1,1-x);
    const double S21mx  = wgplg(2,1,-x);
    const double S111mx = dilog(1-x);
    const double S11mx  = dilog(-x);

    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx2 * lnx;
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx2 * ln1mx;
    const double ln1px  = log(1+x);
    const double ln1px2 = ln1px * ln1px;

    // CF * TR  part
    const double a01 =
      ( 1 - 2 * x + 2 * x2 ) * ( 8 * zeta3 + 4 * ln1mx3 / 3 - 8 * ln1mx * S111mx
				 + 8 * zeta2 * lnx - 4 * lnx * ln1mx2 + 2 * lnx3 / 3 
				 - 8 * lnx * S111mx + 8 * S211mx - 24 * S121mx );
    const double b01 =
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
    const double a02 =
      ( 1 - 2 * x + 2 * x2 ) * ( - 4 * ln1mx3 / 3 + 8 * ln1mx * S111mx - 8 * S211mx ) 
      + ( 1 + 2 * x + 2 * x2 ) * ( - 8 * zeta2 * ln1px - 16 * ln1px * S11mx - 8 * lnx * ln1px2 
				   + 4 * lnx2 * ln1px + 8 * lnx * S11mx - 8 * S21mx - 16 * S12mx )
      + ( 16 + 64 * x ) * ( 2 * S121mx + lnx * S111mx )
      - ( 4 + 8 * x ) * lnx3 / 3 
      + ( 8 - 32 * x + 16 * x2 ) * zeta3 
      - ( 16 + 64 * x ) * zeta2 * lnx;
    const double b02 =
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

  //_________________________________________________________________________________
  ANS2qqH_0::ANS2qqH_0():
    Expression()
  {
  }
  double ANS2qqH_0::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double a0    =
      ( 1 + x2 ) * ( 2 * lnx2 / 3 + 20 * lnx / 9 ) / ( 1 - x ) 
      + 8 * ( 1 - x ) * lnx / 3 + 44 / 27. - 268 * x / 27;
    return CF * TR * a0;
  }
  double ANS2qqH_0::Singular(double const& x) const
  {
    const double a0 = 224 / 27.;
    return CF * TR * a0 / ( 1 - x );
  }
  double ANS2qqH_0::Local(double const& x) const
  {
    const double ln1mx = log(1-x);
    const double a0 = - 8 * zeta3 / 3 + 40 * zeta2 / 9 + 73 / 18. + 224 * ln1mx / 27;
    return CF * TR * a0;
  }

  //_________________________________________________________________________________
  AS2gqH_0::AS2gqH_0():
    Expression()
  {
  }
  double AS2gqH_0::Regular(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double a0    =
      4 * ( 2 / x - 2 + x ) * ln1mx * ln1mx / 3 
      + 8 * ( 10 / x - 10 + 8 * x ) * ln1mx / 9 
      + ( 448 / x - 448 + 344 * x ) / 27;
    return CF * TR * a0;
  }

  //_________________________________________________________________________________
  AS2ggH_0::AS2ggH_0():
    Expression()
  {
  }
  double AS2ggH_0::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double lnx3  = lnx2 * lnx;
    const double ln1mx = log(1-x);
    // CF * TR part
    const double a01   =
      4 * ( 1 + x ) * lnx3 / 3 + ( 6 + 10 * x) * lnx2 
      + ( 32 + 48 * x ) * lnx - 8 / x + 80 - 48 * x - 24 * x2;
    // CA * TR part
    const double a02   =
      4 * ( 1 + x ) * lnx2 / 3 + ( 52 + 88 * x ) * lnx / 9
      - 4 * x * ln1mx / 3 + ( 556 / x - 628 + 548 * x - 700 * x2 ) / 27;
    return TR * ( CF * a01 + CA * a02 );
  }
  double AS2ggH_0::Singular(double const& x) const
  {
    const double a0 = 224 / 27.;
    return CA * TR * a0 / ( 1 - x );
  }
  double AS2ggH_0::Local(double const& x) const
  {
    const double ln1mx = log(1-x);
    const double a01 = - 15.;
    const double a02 = 10 / 9. + 224 * ln1mx / 27;
    return TR * ( CF * a01 + CA * a02 );
  }
}
