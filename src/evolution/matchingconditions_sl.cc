//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingconditions_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"
#include "apfel/integrator.h"

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
  AS1HH_L::AS1HH_L():
    Expression()
  {
  }
  double AS1HH_L::Singular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(x, 2) ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  AS1HH_0::AS1HH_0():
    Expression()
  {
  }
  double AS1HH_0::Singular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(x, 2) ) * ( - 1 - 2 * log(1 - x) ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  AS1gH_L::AS1gH_L():
    Expression()
  {
  }
  double AS1gH_L::Regular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(1 - x, 2) ) / x;
  }

  //_________________________________________________________________________________
  AS1gH_0::AS1gH_0():
    Expression()
  {
  }
  double AS1gH_0::Regular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(1 - x, 2) ) * ( - 1 - 2 * log(x) ) / x;
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
    const double Li21mx = dilog(1 - x);
    const double S121mx = wgplg(1,2,1 - x);
    const double a0 =
      ( 1 + x ) * ( 32 * S121mx + 16 * lnx * Li21mx - 16 * zeta2 * lnx - 4 * pow(lnx,3) / 3 )
      + ( 32 / 3. / x + 8 - 8 * x - 32 * x2 / 3 ) * ( Li21mx - zeta2 )
      + ( 2 + 10 * x + 16 * x2 / 3 ) * lnx * lnx - ( 56 / 3. + 88 * x / 3 + 448 * x2 / 9 ) * lnx
      - 448 / 27. / x - 4 / 3. - 124 * x / 3 + 1600 * x2 / 27;
    return CF * TR * a0;
  }

  //_________________________________________________________________________________
  APS2Hq_L::APS2Hq_L():
    Expression()
  {
  }
  double APS2Hq_L::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double x2    = x * x;
    const double omeL1 =
      8 * ( 1 + x ) * lnx2 - ( 8 + 40 * x + 64 * x2 / 3 ) * lnx
      - 160. / 9 / x + 16 - 48 * x + 448 * x2 / 9;
    return - CF * TR * omeL1;
  }

  //_________________________________________________________________________________
  APS2Hq_L2::APS2Hq_L2():
    Expression()
  {
  }
  double APS2Hq_L2::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double x2    = x * x;
    const double omeL2 = - 8 * ( 1 + x ) * lnx - 16. / 3 / x - 4 + 4 * x + 16 * x2 / 3;
    return CF * TR * omeL2;
  }

  //_________________________________________________________________________________
  AS2Hg_0::AS2Hg_0():
    Expression()
  {
  }
  double AS2Hg_0::Regular(double const& x) const
  {
    const double S121mx = wgplg(1,2,1 - x);
    const double S12mx  = wgplg(1,2,-x);
    const double S211mx = wgplg(2,1,1 - x);
    const double S21mx  = wgplg(2,1,-x);
    const double S111mx = dilog(1 - x);
    const double S11mx  = dilog(-x);

    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx2 * lnx;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx2 * ln1mx;
    const double ln1px  = log(1 + x);
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
  AS2Hg_MSbarMass_0::AS2Hg_MSbarMass_0():
    Expression(),
    _h1(4 * CF),
    _AS2Hg_0(AS2Hg_0{})
  {
  }
  double AS2Hg_MSbarMass_0::Regular(double const& x) const
  {
    return _AS2Hg_0.Regular(x) - 2 * _h1 * 4 * TR * ( pow(x, 2) + pow(1 - x, 2) );
  }

  //_________________________________________________________________________________
  AS2Hg_L::AS2Hg_L():
    Expression()
  {
  }
  double AS2Hg_L::Regular(double const& x) const
  {
    const double x2      = x * x;
    const double lnx     = log(x);
    const double lnx2    = lnx * lnx;
    const double ln1mx   = log(1 - x);
    const double ln1mx2  = ln1mx * ln1mx;
    const double ln1px   = log(1 + x);
    const double S11mx   = dilog(-x);
    const double omeL1p1 =
      CF * TR * ( ( 8 - 16 * x + 16 * x2 ) * ( 2 * lnx * ln1mx - ln1mx2 + 2 * zeta2 )
                  - ( 4 - 8 * x + 16 * x2 ) * lnx2 - 32 * x * ( 1 - x ) * ln1mx
                  - ( 12 - 16 * x + 32 * x2 ) * lnx - 56 + 116 * x - 80 * x2 );
    const double omeL1p2 =
      CA * TR * ( ( 16 + 32 * x + 32 * x2 ) * ( S11mx + lnx * ln1px ) + ( 8 - 16 * x + 16 * x2)
                  * ln1mx2 + ( 8 + 16 * x ) * lnx2 + 32 * x * zeta2 + 32 * x * ( 1 - x ) * ln1mx
                  - ( 8 + 64 * x + 352 * x2 / 3 ) * lnx - 160. / 9 / x + 16 - 200 * x + 1744 * x2 / 9 );
    return - omeL1p1 - omeL1p2;
  }

  //_________________________________________________________________________________
  AS2Hg_MSbarMass_L::AS2Hg_MSbarMass_L():
    Expression(),
    _h1(3 * CF),
    _AS2Hg_L(AS2Hg_L{})
  {
  }
  double AS2Hg_MSbarMass_L::Regular(double const& x) const
  {
    return _AS2Hg_L.Regular(x) - 2 * _h1 * 4 * TR * ( pow(x, 2) + pow(1 - x, 2) );
  }

  //_________________________________________________________________________________
  AS2Hg_L2::AS2Hg_L2():
    Expression()
  {
  }
  double AS2Hg_L2::Regular(double const& x) const
  {
    const double x2      = x * x;
    const double lnx     = log(x);
    const double ln1mx   = log(1 - x);
    const double omeL2p1 =
      CF * TR * ( ( 8 - 16 * x + 16 * x2 ) * ln1mx
                  - ( 4 - 8 * x + 16 * x2 ) * lnx - ( 2 - 8 * x ) );
    const double omeL2p2 =
      CA * TR * ( - ( 8 - 16 * x + 16 * x2 ) * ln1mx - ( 8 + 32 * x ) * lnx
                  - 16. / 3 / x - 4 - 32 * x + 124 * x2 / 3 );
    const double omeL2p3 = TR * TR * ( - 16 * ( x2 + pow(1 - x,2) ) / 3 );
    return omeL2p1 + omeL2p2 + omeL2p3;
  }

  //_________________________________________________________________________________
  ANS2qqH_0::ANS2qqH_0():
    Expression()
  {
  }
  double ANS2qqH_0::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    const double a0   =
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
    const double ln1mx = log(1 - x);
    const double a0    = - 8 * zeta3 / 3 + 40 * zeta2 / 9 + 73 / 18. + 224 * ln1mx / 27;
    return CF * TR * a0;
  }

  //_________________________________________________________________________________
  ANS2qqH_L::ANS2qqH_L():
    Expression()
  {
  }
  double ANS2qqH_L::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double omeL1 = 8 * ( 1 + x2 ) * lnx / 3 / ( 1 - x ) + 8. / 9 - 88 * x / 9;
    return - CF * TR * omeL1;
  }
  double ANS2qqH_L::Singular(double const& x) const
  {
    const double omeL1 = 80. / 9 / ( 1 - x );
    return - CF * TR * omeL1;
  }
  double ANS2qqH_L::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL1 = 80 * ln1mx / 9 + 16 * zeta2 / 3 + 2. / 3;
    return - CF * TR * omeL1;
  }

  //_________________________________________________________________________________
  ANS2qqH_L2::ANS2qqH_L2():
    Expression()
  {
  }
  double ANS2qqH_L2::Regular(double const& x) const
  {
    const double omeL2 = - 4. / 3 - 4 * x / 3;
    return CF * TR * omeL2;
  }
  double ANS2qqH_L2::Singular(double const& x) const
  {
    const double omeL2 = 8. / 3. / ( 1 - x );
    return CF * TR * omeL2;
  }
  double ANS2qqH_L2::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL2 = 8 * ln1mx / 3 + 2;
    return CF * TR * omeL2;
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
  AS2gqH_L::AS2gqH_L():
    Expression()
  {
  }
  double AS2gqH_L::Regular(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL1 =
      160. / 9 / x - 160. / 9 + 128 * x / 9
      + ( 32. / 3 / x - 32. / 3 + 16 * x / 3 ) * ln1mx;
    return - CF * TR * omeL1;
  }

  //_________________________________________________________________________________
  AS2gqH_L2::AS2gqH_L2():
    Expression()
  {
  }
  double AS2gqH_L2::Regular(double const& x) const
  {
    const double omeL2 = 16. / 3 / x - 16. / 3 + 8 * x / 3;
    return CF * TR * omeL2;
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
    const double ln1mx = log(1 - x);
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
    const double ln1mx = log(1 - x);
    const double a01   = - 15.;
    const double a02   = 10 / 9. + 224 * ln1mx / 27;
    return TR * ( CF * a01 + CA * a02 );
  }

  //_________________________________________________________________________________
  AS2ggH_MSbarMass_0::AS2ggH_MSbarMass_0():
    Expression(),
    _h1(4 * CF),
    _AS2ggH_0(AS2ggH_0{})
  {
  }
  double AS2ggH_MSbarMass_0::Regular(double const& x) const
  {
    return _AS2ggH_0.Regular(x);
  }
  double AS2ggH_MSbarMass_0::Singular(double const& x) const
  {
    return _AS2ggH_0.Singular(x);
  }
  double AS2ggH_MSbarMass_0::Local(double const& x) const
  {
    return _AS2ggH_0.Local(x) - 2 * _h1 * ( - 4. * TR / 3. );
  }

  //_________________________________________________________________________________
  AS2ggH_L::AS2ggH_L():
    Expression()
  {
  }
  double AS2ggH_L::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double omeL1 =
      + CF * TR * ( 8 * ( 1 + x ) * lnx2 + ( 24 + 40 * x ) * lnx
                    - 16. / 3 / x + 64 - 32 * x - 80 * x2 / 3 )
      + CA * TR * ( 16 * ( 1 + x ) * lnx / 3 + 184. / 9 / x
                    - 232. / 9 + 152 * x / 9 - 184 * x2 / 9 );
    return - omeL1;
  }
  double AS2ggH_L::Singular(double const& x) const
  {
    const double omeL1 = 80. / 9 / ( 1 - x );
    return - CA * TR * omeL1;
  }
  double AS2ggH_L::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL1 = CF * TR * 4 + CA * TR * ( 16. / 3 + 80 * ln1mx / 9 );
    return - omeL1;
  }

  //_________________________________________________________________________________
  AS2ggH_MSbarMass_L::AS2ggH_MSbarMass_L():
    Expression(),
    _h1(3 * CF),
    _AS2ggH_L(AS2ggH_L{})
  {
  }
  double AS2ggH_MSbarMass_L::Regular(double const& x) const
  {
    return _AS2ggH_L.Regular(x);
  }
  double AS2ggH_MSbarMass_L::Singular(double const& x) const
  {
    return _AS2ggH_L.Singular(x);
  }
  double AS2ggH_MSbarMass_L::Local(double const& x) const
  {
    return _AS2ggH_L.Local(x) - 2 * _h1 * ( - 4. * TR / 3. );
  }

  //_________________________________________________________________________________
  AS2ggH_L2::AS2ggH_L2():
    Expression()
  {
  }
  double AS2ggH_L2::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double omeL2 =
      + CF * TR * ( 8 * ( 1 + x ) * lnx + 16. / 3 / x + 4
                    - 4 * x - 16 * x2 / 3 )
      + CA * TR * ( 8. / 3 / x - 16. / 3 + 8 * x / 3
                    - 8 * x2 / 3 );
    return omeL2;
  }
  double AS2ggH_L2::Singular(double const& x) const
  {
    const double omeL2 = 8. / 3 / ( 1 - x );
    return CA * TR * omeL2;
  }
  double AS2ggH_L2::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    const double omeL2 = TR * TR * 16. / 9 + CA * TR * 8 * ln1mx / 3;
    return omeL2;
  }

  //_________________________________________________________________________________
  APS3Hq_0::APS3Hq_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double APS3Hq_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return ps_(&xr, &nfr, &asr, &LLr) + (x < 0.5 ? aps1_(&xr, &nfr, &asr, &LLr) : aps2_(&xr, &nfr, &asr, &LLr));
  }

  //_________________________________________________________________________________
  AS3Hg_0::AS3Hg_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double AS3Hg_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return qg_(&xr, &nfr, &asr, &LLr) + aqg3_(&xr) / 2 + (x < 0.5 ? red0_(&xr) : red1_(&xr)) + _nf * aqg3nf_(&xr);
  }

  //_________________________________________________________________________________
  ANS3qqH_0::ANS3qqH_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double ANS3qqH_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return nsreg_(&xr, &nfr, &asr, &LLr) + ansreg_(&xr, &nfr, &asr, &LLr);
  }
  double ANS3qqH_0::Singular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return nsplu_(&xr, &nfr, &asr, &LLr) + ansplu1_(&xr, &nfr, &asr, &LLr) + ansplu2_(&xr, &nfr, &asr, &LLr);
  }
  double ANS3qqH_0::Local(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    const double Iansplu1 = Integrator{[=] (double const& y) -> double
      {
        double yr  = y;
        double nfr = _nf;
        double asr = 1;
        double LLr = 0;
        return ansplu1_(&yr, &nfr, &asr, &LLr);
      }
    }.integrate(0, x, eps5);
    const double Iansplu2 = Integrator{[=] (double const& y) -> double
      {
        double yr  = y;
        double nfr = _nf;
        double asr = 1;
        double LLr = 0;
        return ansplu2_(&yr, &nfr, &asr, &LLr);
      }
    }.integrate(0, x, eps5);
    return nsdel_(&xr, &nfr, &asr, &LLr) + ansdel_(&xr, &nfr, &asr, &LLr) + nsplu_(&xr, &nfr, &asr, &LLr) * ( 1 - x ) * log(1 - x) - Iansplu1 - Iansplu2;
  }

  //_________________________________________________________________________________
  AS3gqH_0::AS3gqH_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double AS3gqH_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return gq_(&xr, &nfr, &asr, &LLr) + agq_(&xr, &nfr, &asr, &LLr);
  }

  //_________________________________________________________________________________
  AS3ggH_0::AS3ggH_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double AS3ggH_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return ggreg_(&xr, &nfr, &asr, &LLr) + (x < 0.5 ? aggreg0_(&xr, &nfr, &asr, &LLr) : aggreg1_(&xr, &nfr, &asr, &LLr));
  }
  double AS3ggH_0::Singular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return ggplu_(&xr, &nfr, &asr, &LLr) + aggplu_(&xr, &nfr, &asr, &LLr);
  }
  double AS3ggH_0::Local(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return ggdel_(&xr, &nfr, &asr, &LLr) + aggdel_(&xr, &nfr, &asr, &LLr) + ( ggplu_(&xr, &nfr, &asr, &LLr) + aggplu_(&xr, &nfr, &asr, &LLr) ) * ( 1 - x ) * log(1 - x);
  }

  //_________________________________________________________________________________
  AS3qgQ_0::AS3qgQ_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double AS3qgQ_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return qgl_(&xr, &nfr, &asr, &LLr);
  }

  //_________________________________________________________________________________
  APS3qqQ_0::APS3qqQ_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double APS3qqQ_0::Regular(double const& x) const
  {
    double xr  = x;
    double nfr = _nf;
    double asr = 1;
    double LLr = 0;
    return psl_(&xr, &nfr, &asr, &LLr);
  }

  //_________________________________________________________________________________
  AS1polHg_L::AS1polHg_L():
    Expression()
  {
  }
  double AS1polHg_L::Regular(double const& x) const
  {
    return 4 * TR * ( 2 * x - 1 );
  }

  //_________________________________________________________________________________
  AS1polggH_L::AS1polggH_L():
    Expression()
  {
  }
  double AS1polggH_L::Local(double const&) const
  {
    return - 4 * TR / 3.;
  }

  //_________________________________________________________________________________
  AS1polHH_L::AS1polHH_L():
    Expression()
  {
  }
  double AS1polHH_L::Singular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(x, 2) ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  AS1polHH_0::AS1polHH_0():
    Expression()
  {
  }
  double AS1polHH_0::Singular(double const& x) const
  {
    return 2 * CF * ( 1 + pow(x, 2) ) * ( - 1 - 2 * log(1 - x) ) / ( 1 - x );
  }

  //_________________________________________________________________________________
  AS1polgH_L::AS1polgH_L():
    Expression()
  {
  }
  double AS1polgH_L::Regular(double const& x) const
  {
    return 2 * CF * ( 2 - x );
  }

  //_________________________________________________________________________________
  AS1polgH_0::AS1polgH_0():
    Expression()
  {
  }
  double AS1polgH_0::Regular(double const& x) const
  {
    return 2 * CF * ( 2 - 2 * x - 2 * ( 2 - x ) * log(x) );
  }

  //_________________________________________________________________________________
  APS2polHq_0::APS2polHq_0(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double APS2polHq_0::Regular(double const& x) const
  {
    // Eq. (A.4) of https://arxiv.org/pdf/hep-ph/9608342
    const double S121mx = wgplg(1, 2, 1 - x);
    const double Li21mx = dilog(1 - x);
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    return CF * TR * ( ( 1 + x ) * ( 32 * S121mx + 16 * lnx * Li21mx - 24 * zeta2 * lnx - 4 * lnx3 / 3 )
                       + 20 * ( 1 - x ) * ( 2 * Li21mx - 3 * zeta2 ) - ( 2 - 6 * x ) * lnx2
                       - ( 12 + 60 * x ) * lnx - 72 * ( 1 - x ) )
           // Last term in Eq. (138) of https://arxiv.org/pdf/2211.15337
           + 4 * CF * TR * zeta2 * ( 5 * ( 1 - x ) + 2 * ( 1 + x ) * lnx )
           // Second term in the r.h.s. of Eq. (141) of https://arxiv.org/pdf/2211.15337
           + 4 * CF * TR * _nf * ( ( 2 + x ) * lnx2 + 2 * ( 3 - x ) * lnx + 4 * ( 1 - x ) );
  }

  //_________________________________________________________________________________
  APS2polHq_L::APS2polHq_L():
    Expression()
  {
  }
  double APS2polHq_L::Regular(double const& x) const
  {
    const double lnx = log(x);
    return - 8 * CF * TR * ( - 1 + x + lnx * ( 1 - 3 * x + ( 1 + x ) * lnx ) );
  }

  //_________________________________________________________________________________
  APS2polHq_L2::APS2polHq_L2():
    Expression()
  {
  }
  double APS2polHq_L2::Regular(double const& x) const
  {
    return - 4 * CF * TR * ( 5 - 5 * x + 2 * ( 1 + x ) * log(x) );
  }

  //_________________________________________________________________________________
  AS2polHg_0::AS2polHg_0():
    Expression()
  {
  }
  double AS2polHg_0::Regular(double const& x) const
  {
    // Eq. (A.2) of https://arxiv.org/pdf/hep-ph/9608342
    const double x2     = x * x;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double ln1px  = log(1 + x);
    const double ln1px2 = ln1px * ln1px;
    const double Li21mx = dilog(1 - x);
    const double Li2mx  = dilog(- x);
    const double Li31mx = trilog(1 - x);
    const double Li3mx  = trilog(- x);
    const double S121mx = wgplg(1, 2, 1 - x);
    const double S12mx  = wgplg(1, 2, - x);
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    return - ( CF * TR * ( ( - 1 + 2 * x ) * ( 8 * zeta3 + 8 * zeta2 * ln1mx + 4 * ln1mx3 / 3 - 8 * ln1mx * Li21mx + 4 * zeta2 * lnx - 4 * lnx * ln1mx2
                                               + 2 * lnx3 / 3 - 8 * lnx * Li21mx + 8 * Li31mx - 24 * S121mx )
                           - ( 116 - 48 * x - 16 * x2 ) * Li21mx + ( 50 - 32 * x - 8 * x2 ) * zeta2 - ( 72 - 16 * x - 8 * x2 ) * lnx * ln1mx
                           + ( 12 - 8 * x - 4 * x2 ) * ln1mx2 - ( 5 - 8 * x - 4 * x2 ) * lnx2 - ( 64 - 60 * x ) * ln1mx - ( 16 + 50 * x ) * lnx - 22 + 46 * x )
               + CA * TR * ( ( - 1 + 2 * x ) * ( - 8 * zeta2 * ln1mx - 4 * ln1mx3 / 3 + 8 * ln1mx * Li21mx - 8 * Li31mx )
                             + ( 1 + 2 * x) * ( - 4 * lnx3 / 3 - 8 * zeta2 * ln1px - 16 * ln1px * Li2mx - 8 * lnx * ln1px2 + 4 * lnx2 * ln1px + 8 * lnx * Li2mx
                                                - 8 * Li3mx - 16 * S12mx )
                             + 16 * ( 1 + x ) * ( 4 * S121mx + 2 * lnx * Li21mx - 3 * zeta2 * lnx + Li2mx + lnx * ln1px ) - 16 * ( 1 - x ) * zeta3
                             + ( 100 - 112 * x - 8 * x2 ) * Li21mx - ( 132 - 144 * x - 4 * x2 ) * zeta2 - 4 * x * ( 4 +  x ) * lnx * ln1mx
                             - ( 10 - 8 * x - 2 * x2 ) * ln1mx2 - ( 6 + 2 * x2 ) * lnx2 + 4 * ln1mx - ( 56 + 148 * x ) * lnx - 204 + 212 * x ) )
           // Last term in the second line of Eq. (111) of https://arxiv.org/pdf/2211.15337
           + 2 * CF * TR * ( - 3 - 2 * ( 1 - 2 * x ) * lnx + 4 * ( 1 - 2 * x ) * ln1mx )
           + 8 * CA * TR * ( 6 * ( 1 - x ) + 2 * ( 1 + x ) * lnx - ( 1 - 2 * x ) * ln1mx );
  }

  //_________________________________________________________________________________
  AS2polHg_L::AS2polHg_L():
    Expression()
  {
  }
  double AS2polHg_L::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double ln1mx = log(1 - x);
    const double ln1px = log(1 + x);
    return - 4. / 3 * TR * ( 8 * TR * ( 1 - 2 * x ) + 66 * CA * x + 12 * CA * ( - 6 + zeta2 )
                             + 3 * CF * ( - 2 - 4 * zeta2 + x * ( - 3 + 8 * zeta2 ) )
                             + 6 * ( CA - CF ) * ( - 1 + 2 * x ) * pow(ln1mx, 2)
                             + 3 * ln1mx * ( - 8 * ( CA - CF ) * ( - 1 + x ) + 4 * CF * ( - 1 + 2 * x ) * lnx )
                             + 3 * lnx * ( CF - 2 * ( CA + 8 * ( CA + CF ) * x ) + ( CF - 2 * CF * x + CA * ( 2 + 4 * x ) ) * lnx
                                           + 4 * ( CA + 2 * CA * x ) *ln1px ) + 12 * CA * ( 1 + 2 * x ) * dilog(-x) );
  }

  //_________________________________________________________________________________
  AS2polHg_L2::AS2polHg_L2():
    Expression()
  {
  }
  double AS2polHg_L2::Regular(double const& x) const
  {
    return - 2 * TR * ( - 3 * ( CF + 8 * CA * ( - 1 + x ) ) + 4 * ( CA - CF ) * ( - 1 + 2 * x ) * log(1 - x)
                        + ( - 2 * CF + 4 * CF * x + 8 * CA * ( 1 + x ) ) * log(x) );
  }

  //_________________________________________________________________________________
  ANS2polqqH_0::ANS2polqqH_0():
    Expression()
  {
  }
  double ANS2polqqH_0::Regular(double const& x) const
  {
    return 2. / 27 * CF * TR * ( 22 - 134 * x - 3 * log(x) * ( 22 - 24 * x + 22 * pow(x, 2) + 3 * ( 1 + pow(x, 2) ) * log(x) ) / ( - 1 + x ) );
  }
  double ANS2polqqH_0::Singular(double const& x) const
  {
    return CF * TR * 224. / 27 / ( 1 - x );
  }
  double ANS2polqqH_0::Local(double const& x) const
  {
    return CF * TR / 18 * ( 73 + 80 * zeta2 - 48 * zeta3 ) + CF * TR * 224. / 27 * log(1 - x);
  }

  //_________________________________________________________________________________
  ANS2polqqH_L::ANS2polqqH_L():
    Expression()
  {
  }
  double ANS2polqqH_L::Regular(double const& x) const
  {
    return 8. / 9 * CF * TR * ( - 1 + 11 * x + 3 * ( 1 + pow(x, 2) ) * log(x) / ( - 1 + x ) );
  }
  double ANS2polqqH_L::Singular(double const& x) const
  {
    return - CF * TR * 80. / 9 / ( 1 - x );
  }
  double ANS2polqqH_L::Local(double const& x) const
  {
    return - CF * TR / 18 *( 12 + 96 * zeta2 ) - CF * TR * 80. / 9 * log(1 - x);
  }

  //_________________________________________________________________________________
  ANS2polqqH_L2::ANS2polqqH_L2():
    Expression()
  {
  }
  double ANS2polqqH_L2::Regular(double const& x) const
  {
    return -4. / 3 * CF * TR * ( 1 + x );
  }
  double ANS2polqqH_L2::Singular(double const& x) const
  {
    return CF * TR * 8. / 3 / ( 1 - x );
  }
  double ANS2polqqH_L2::Local(double const& x) const
  {
    return 2 * CF * TR + CF * TR * 8. / 3 * log(1 - x);
  }

  //_________________________________________________________________________________
  AS2polgqH_0::AS2polgqH_0():
    Expression()
  {
  }
  double AS2polgqH_0::Regular(double const& x) const
  {
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    return CF * TR * ( - ( 32 * ( - 11 + 4 * x ) ) / 27 + ( 8 * ( 4 + x ) * ln1mx ) / 9 - ( 4 * ( - 2 + x ) * ln1mx2 ) / 3 );
  }

  //_________________________________________________________________________________
  AS2polgqH_L::AS2polgqH_L():
    Expression()
  {
  }
  double AS2polgqH_L::Regular(double const& x) const
  {
    return - CF * TR * ( 16 * ( 4 + x ) / 9 - 16 * ( - 2 + x ) * log(1 - x) / 3 );
  }

  //_________________________________________________________________________________
  AS2polgqH_L2::AS2polgqH_L2():
    Expression()
  {
  }
  double AS2polgqH_L2::Regular(double const& x) const
  {
    return - CF * TR * 8 * ( - 2 + x ) / 3;
  }

  //_________________________________________________________________________________
  AS2polggH_0::AS2polggH_0():
    Expression()
  {
  }
  double AS2polggH_0::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double lnx3  = lnx * lnx2;
    const double ln1mx = log(1 - x);
    return CA * TR * ( - 2 * ( - 337 + 449 * x ) / 27 - 4 * x * ln1mx / 3 + 4 * ( 22 + x ) * lnx / 9  + 4 * ( 1 + x ) * lnx2 / 3 )
           + CF * TR * ( - 56 * ( - 1 + x ) + 12 * ( 3 + x ) * lnx -  2 * ( - 5 + x ) * lnx2 + 4 * ( 1 + x ) * lnx3 / 3 );
  }
  double AS2polggH_0::Singular(double const& x ) const
  {
    return CA * 224. / 27 * TR / ( 1 - x );
  }
  double AS2polggH_0::Local(double const& x) const
  {
    return - CF * 15 * TR + 2. / 9 * CA * 5 * TR + CA * 224. / 27 * TR * log(1 - x);
  }

  //_________________________________________________________________________________
  AS2polggH_L::AS2polggH_L():
    Expression()
  {
  }
  double AS2polggH_L::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return - CA * TR * ( - 16 * ( - 14 + 19 * x ) / 9 + 16 * ( 1 + x ) * lnx / 3 )
           - CF * TR * ( - 40 * ( - 1 + x ) - 8 * ( - 5 + x ) * lnx + 8 * ( 1 + x ) * lnx2 );
  }
  double AS2polggH_L::Singular(double const& x) const
  {
    return - CA * 80. / 9 * TR / ( 1 - x );
  }
  double AS2polggH_L::Local(double const& x) const
  {
    return - CF * 4 * TR - 2. / 9 * CA * 24 * TR - CA * 80. / 9 * TR * log(1 - x);
  }

  //_________________________________________________________________________________
  AS2polggH_L2::AS2polggH_L2():
    Expression()
  {
  }
  double AS2polggH_L2::Regular(double const& x) const
  {
    const double lnx = log(x);
    return - CA * TR * 8 * ( - 1 + 2 * x ) / 3 + CF * TR * ( - 20 * ( - 1 + x ) + 8 * ( 1 + x ) * lnx );
  }
  double AS2polggH_L2::Singular(double const& x) const
  {
    return CA * 8. / 3 * TR / ( 1 - x );
  }
  double AS2polggH_L2::Local(double const& x) const
  {
    return 16 * pow(TR, 2) / 9 + CA * 8. / 3 * TR * log(1 - x);
  }
}
