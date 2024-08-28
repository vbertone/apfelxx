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
  AS2Hg_L::AS2Hg_L():
    Expression()
  {
  }
  double AS2Hg_L::Regular(double const& x) const
  {
    const double x2      = x * x;
    const double lnx     = log(x);
    const double lnx2    = lnx * lnx;
    const double ln1mx   = log(1-x);
    const double ln1mx2  = ln1mx * ln1mx;
    const double ln1px   = log(1+x);
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
  AS2Hg_L2::AS2Hg_L2():
    Expression()
  {
  }
  double AS2Hg_L2::Regular(double const& x) const
  {
    const double x2      = x * x;
    const double lnx     = log(x);
    const double ln1mx   = log(1-x);
    const double omeL2p1 =
      CF * TR * ( ( 8 - 16 * x + 16 * x2 ) * ln1mx
                  - ( 4 - 8 * x + 16 * x2 ) * lnx - ( 2 - 8 * x ) );
    const double omeL2p2 =
      CA * TR * ( - ( 8 - 16 * x + 16 * x2 ) * ln1mx - ( 8 + 32 * x ) * lnx
                  - 16. / 3 / x - 4 - 32 * x + 124 * x2 / 3 );
    const double omeL2p3 = TR * TR * ( - 16 * ( x2 + pow(1-x,2) ) / 3 );
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
    const double ln1mx = log(1-x);
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
    const double ln1mx = log(1-x);
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
    const double ln1mx = log(1-x);
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
    const double ln1mx = log(1-x);
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
    const double a01   = - 15.;
    const double a02   = 10 / 9. + 224 * ln1mx / 27;
    return TR * ( CF * a01 + CA * a02 );
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
    const double ln1mx = log(1-x);
    const double omeL1 = CF * TR * 4 + CA * TR * ( 16. / 3 + 80 * ln1mx / 9 );
    return - omeL1;
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
    const double ln1mx = log(1-x);
    const double omeL2 = TR * TR * 16. / 9 + CA * TR * 8 * ln1mx / 3;
    return omeL2;
  }

  //_________________________________________________________________________________
  APS3Hq_0::APS3Hq_0():
    Expression()
  {
  }
  double APS3Hq_0::Regular(double const& x) const
  {
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double lnx3   = lnx * lnx2;
    const double lnx4   = lnx * lnx3;
    const double lnx5   = lnx * lnx4;
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double ln1mx4 = ln1mx * ln1mx3;
    return pow(1 - x, 2) * ( - 158.3709570235676 * ln1mx3 + 2151.0528216721095 * ln1mx2 - 14178.20613565087 )
           + 0.7136042069623996 * x + 15618.788026509601 * x * lnx2 + 2117.261727957378 * x * lnx + 688.396 * lnx / x
           + ( 1 - x ) * 3812.8990 / x + 1.6 * lnx5 - 20.345678 * lnx4 + 165.11455 * lnx3 - 604.63554 * lnx2 + 3524.9967 * lnx
           + ( 1 - x ) * ( 0.24691358 * ln1mx4 - 4.44444444 * ln1mx3 - 2.28230742 * ln1mx2 - 357.426943 * ln1mx + 116.478169 );
  }

  //_________________________________________________________________________________
  AS3Hg_0::AS3Hg_0(double const& rho):
    Expression(),
    _rho(rho)
  {
    // Moments for the known exact small-x contribution (Vogt)
    const std::vector<double> N{1996.50147366 - _rho, 232.27499707 - _rho / 3, 7.44963895 - _rho / 5,
                                -77.93698013 - _rho * 0.14285714285714285, -121.36688921 - _rho / 9};

    // Matrix
    const std::vector<std::vector<double>>  inv_A
    {
      {-1.13879194e-01, -7.30104359e-01, -1.16019755e+02, 1.47898113e+02, -2.48866805e+01},
      {2.60826752e+00, 1.65606748e+01, 1.87421683e+03, -2.53062040e+03, 2.76071734e+02},
      {-1.16231560e+01, -7.31022502e+01, -6.56760263e+03, 9.21295249e+03, -8.50419108e+02},
      {1.73021125e+01, 1.07828576e+02, 8.21779251e+03, -1.18144214e+04,  9.95444799e+02},
      {-8.20269840e+00, -5.06719525e+01, -3.41032145e+03, 4.98742337e+03, -3.96388523e+02}
    };

    // Matrix multiplication of inv_A and N
    _C.resize(N.size(), 0.);
    for (int i = 0; i < (int) N.size(); i++)
      for (int j = 0; j < (int) N.size(); j++)
        _C[i] += inv_A[j][i] * N[j];
  }
  double AS3Hg_0::Regular(double const& x) const
  {
    if (x > 0.924)
      return Regular(0.924 - eps5) * ( 1 - x );
    else
      return _C[0] * pow(log(1-x), 5) + _C[1] * pow(log(1-x), 4) + _C[2] * x + _C[3] * pow(x, 2) + _C[4] * log(x) + _rho / x
             + ( 41984. / 27. + 160 * zeta2 - 224 * zeta3 ) * log(x) / x;
  }

  //_________________________________________________________________________________
  ANS3qqH_0::ANS3qqH_0(double const& rho):
    Expression(),
    _rho(rho)
  {
    // Moments for the known exact small-x contribution (Vogt)
    const std::vector<double> N{-31.101 + _rho * 0.375, -27.989 + _rho * 0.023437499998337614, -20.146 + _rho * 0.004629629629528962,
                                -13.313 + _rho * 0.0014648437499999065, -7.629 + _rho * 0.0005999999999999996, -2.854 + _rho * 0.0002893518518518518,
                                1.228 + _rho * 0.00015618492294877136};

    // Matrix
    const std::vector<std::vector<double>>  inv_A
    {
      {2.16367103e+01, 2.20120178e+02, -7.28098793e-01, 1.08859149e+03, 3.16215872e+01, 7.57356352e+02, -1.31967754e+02},
      {-1.30805003e+03, -1.30363921e+04, 4.71551230e+01, -6.03760425e+04, -1.06557895e+03, -4.38642686e+04, 5.99320146e+03},
      {1.36620524e+04,  1.34037998e+05, -5.20816646e+02, 5.95121449e+05, 8.47026017e+03, 4.43696299e+05, -5.19148509e+04},
      {-5.23812588e+04, -5.07554870e+05,2.09319457e+03, -2.18873394e+06, -2.74098723e+04, -1.65937417e+06, 1.75138916e+05},
      {9.13767640e+04, 8.76404430e+05, -3.80468317e+03, 3.69828317e+06, 4.24347906e+04, 2.83722814e+06, -2.77933570e+05},
      {-7.40399878e+04, -7.04036144e+05, 3.19808122e+03,-2.92103432e+06, -3.14155811e+04, -2.26097908e+06,2.09168184e+05},
      {2.26683150e+04, 2.13961524e+05, -1.01233621e+03, 8.75640617e+05, 8.95431572e+03, 6.82526398e+05, -6.03194325e+04}
    };

    // Matrix multiplication of inv_A and N
    _C.resize(N.size(), 0.);
    for (int i = 0; i < (int) N.size(); i++)
      for (int j = 0; j < (int) N.size(); j++)
        _C[i] += inv_A[j][i] * N[j];
  }
  double ANS3qqH_0::Regular(double const& x) const
  {
    if (x > 0.4)
      return 0;
    else
      return _C[0] * pow(log(1-x), 3) + _C[1] * pow(log(1-x), 2) +_C[3] * x + _C[4] * pow(log(x), 2) + _C[5] * log(1-x) + _C[6] + _rho * pow(log(x), 3);
  }
  double ANS3qqH_0::Singular(double const& x) const
  {
    if (x > 0.4)
      return 0;
    else
      return _C[2] / ( 1 - x );
  }
  double ANS3qqH_0::Local(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  AS3gqH_0::AS3gqH_0():
    Expression()
  {
  }
  double AS3gqH_0::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1-x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double ln1mx4 = ln1mx * ln1mx3;
    if (x > 0.7)
      return Regular(0.7 - eps5) * ( 1 - x ) / 0.3;
    else
      return - 237.1720947621626 * ln1mx3 - 201.496990873891 * ln1mx2 + 7247.6979315802555 * ln1mx + 39967.310031716734 * x2
             - 22017.709314369586 - 28459.052011351927 * lnx - 14511.477130954034 * lnx2 - 341.543 * lnx / x
             + 1814.73 / x - 580. / 243. * ln1mx4 - 17624. / 729. * ln1mx3 - 135.699 * ln1mx2;
  }

  //_________________________________________________________________________________
  AS3ggH_0::AS3ggH_0(double const& rho):
    Expression(),
    _rho(rho)
  {
    // Moments for the known exact small-x contribution (Vogt)
    const std::vector<double> N{-441.3444 + _rho * 0.9999999999999999, -96.32 + _rho * 0.11111111111076043,
                                -13.33 + _rho * 0.040000000000137814, 31.1 + _rho * 0.020408163265306242, 60.84 + _rho * 0.01234567901234568};

    // Matrix
    const std::vector<std::vector<double>>  inv_A
    {
      {1.43759679e+01, -3.46333748e+02, 1.61699126e+03, -2.51174722e+03, 1.23814937e+03},
      {1.21751771e+02, -2.87863393e+03, 1.32164568e+04, -2.02263161e+04, 9.83890667e+03},
      {2.91559485e+02, -5.88750691e+03, 2.44589984e+04, -3.49174315e+04, 1.61247555e+04},
      {-2.80668008e+01, 3.50042899e+02, -1.18494800e+03, 1.50037935e+03, -6.38918620e+02},
      {-3.82520602e+01, 1.94881955e+01, 2.01486391e+03, -5.00918023e+03, 3.06450021e+03}
    };

    // Matrix multiplication of inv_A and N
    _C.resize(N.size(), 0.);
    for (int i = 0; i < (int) N.size(); i++)
      for (int j = 0; j < (int) N.size(); j++)
        _C[i] += inv_A[i][j] * N[j];
  }
  double AS3ggH_0::Regular(double const& x) const
  {
    return _C[0] * pow(log(1-x), 2) + _C[1] * log(1-x) + _C[2] * pow(x, 2) + _C[3] * log(x) + _C[4] * x + _rho * log(x) / x;
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
  APS2polHq_0::APS2polHq_0():
    Expression()
  {
  }
  double APS2polHq_0::Regular(double const& x) const
  {
    const double lnx = log(x);
    return 2. / 3 * CF * TR * ( 84 * ( - 1 + x ) + 48 * ( 1 + x ) * zeta3
                                + lnx * ( 18 - 102 * x + 60 * ( - 1 + x ) * log(1 - x) + lnx * ( 9 + 15 * x - 2 * ( 1 + x ) * lnx ) )
                                + 12 * ( 5 *( - 1 + x ) + 2 * ( 1 + x ) * lnx ) * dilog(x) - 48 * ( 1 + x ) * trilog(x) );
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
    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1mx  = log(1 - x);
    const double ln1mx2 = ln1mx * ln1mx;
    const double ln1mx3 = ln1mx * ln1mx2;
    const double ln1px  = log(1 + x);
    const double ln1px2 = ln1px * ln1px;
    const double atan   = ( log(1 - x) - lnx ) / 2;
    const double S12mx  = wgplg(1,2,-x);
    return TR / 3 * ( 6 * ( CF * ( - 11 - 36 * zeta2 + x * ( 23 + 4 * ( 2 + x ) * zeta2 - 16 * zeta3 ) + 8 * zeta3 )
                            + 2 * CA * ( - 51 + 4 * zeta2 + 12 * zeta3 + x * ( 53 - ( 4 + x ) * zeta2 + 20 * zeta3 ) ) ) - 4 * ( CA - CF ) * ( - 1 + 2 * x ) * ln1mx3
                      + 6 * ln1mx2 * ( ( - 1 + x ) * ( - 2 * CF * ( 3 + x ) + CA * ( 5 + x ) ) - 2 * ( 2 * CA - CF ) * ( - 1 + 2 * x ) * lnx )
                      + lnx * ( - 6 * ( CF * ( 8 + 25 * x ) + CA * ( 28 + 74 * x ) ) - 3 * ( CA * ( 6 + 4 * x ) + CF * ( 5 - 4 * x * ( 2 + x ) ) ) * lnx
                                - 2 * ( CF - 2 * CF * x + CA * ( 2 + 4 * x ) ) * lnx2 )
                      + 4 * ln1mx * ( CF * ( - 48 + Pi2 * ( 1 - 2 * x ) + 45 * x ) + CA * ( 3 + Pi2 * ( - 1 + 2 * x ) )
                                      + 3 * lnx * ( CA * ( - 1 + x ) * ( 25 + x ) + CF * ( 11 - 2 * x * ( 4 + x ) ) + ( CF - 2 * CF * x ) * lnx ) )
                      + 12 * CA * ( - 2 * ( zeta2 + 2 * x * zeta2 ) + lnx * ( 4 + 4 * x + lnx + 2 * x * lnx ) ) * ln1px - 24 * CA * ( 1 + 2 * x ) * lnx * ln1px2
                      + 6 * ( 2 * ( CF * ( 29 - 4 * x * ( 3 + x ) ) + CA * ( - 29 + 2 * x * ( 12 + x ) ) - 8 * CA * x * atan + 2 * ( CA + CF * ( - 1 + 2 * x ) ) * ln1mx
                                    + 2 * ( 3 * CA + 2 * CF - 4 * CF * x ) * lnx + 4 * CA * ( 1 + 2 * x ) * ln1px ) * dilog(x)
                              + 2 * CA * ( 2 + 2 * x + lnx + 2 * x * lnx - 2 * ( 1 + 2 * x ) * ln1px ) * dilog(x2)
                              - 4 * ( CA - CF ) * ( - 1 + 2 * x ) * trilog(1 - x) - 4 * ( CF * ( 3 - 6 * x ) + CA * ( 7 + 6 * x ) ) * trilog(x)
                              - ( CA + 2 * CA * x ) * trilog(x2) ) - 48 * ( CA + 2 * CA * x ) * S12mx );
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
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
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
