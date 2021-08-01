//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingfunctionsff.h"
#include "apfel/constants.h"

#include <numeric>

namespace apfel
{
  //_________________________________________________________________________________
  C1nsff::C1nsff():
    Expression()
  {
  }
  double C1nsff::Regular(double const& x) const
  {
    return CF * ( 2 * ( 1 - x ) + 4 * ( 1 + x * x ) * log(x) / ( 1 - x ) );
  }
  double C1nsff::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  C1qgff::C1qgff():
    Expression()
  {
  }
  double C1qgff::Regular(double const& x) const
  {
    return 2 * CF * ( 2 * x + 4 * ( 1 + pow(1 - x, 2) ) * log(x) / x );
  }

  //_________________________________________________________________________________
  C1gqff::C1gqff():
    Expression()
  {
  }
  double C1gqff::Regular(double const& x) const
  {
    return TR * ( 4 * x * ( 1 - x ) + 4 * ( 1 - 2 * x * ( 1 - x ) ) * log(x) );
  }

  //_________________________________________________________________________________
  C1ggff::C1ggff():
    Expression()
  {
  }
  double C1ggff::Regular(double const& x) const
  {
    return 8 * CA * pow(1 - x * ( 1 - x ), 2) * log(x) / x / ( 1 - x );
  }
  double C1ggff::Local(double const&) const
  {
    return - CA * zeta2;
  }

  //_________________________________________________________________________________
  C2nspff::C2nspff(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = 14.926669450170849 + 5.530864197530864*_nf;
  }
  double C2nspff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double x4   = x * x3;
    const double x5   = x * x4;
    const double x6   = x * x5;
    const double xb   = 1 - x;
    const double xb2  = xb * xb;
    const double xb3  = xb * xb2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    return -373.6403274234052 - 190.23261022266544*x4 + 11.62459533324049*x5 - 0.48152936575476524*x6
           - 22.22222222222222*l1x + 7.111111111111111*l1x2
           + xb*(141.04140553909332 + 47.111111111111114*l1x - 3.5555555555555554*l1x2)
           + xb3*(3.9346330275229358*l1x - 0.7608247422680412*l1x2)
           + xb2*(-12.865738999623918*l1x + 2.426546391752577*l1x2)
           - 40.5736752705563*lx - 32.66666666666667*lx2 + 3.259259259259259*lx3
           + (2*(-1.4599875154038369 + 3.5555555555555554*lx + 7.111111111111111*lx2))/x
           + x2*(56.065623348965666 - 30.13516702067784*lx - 7.172955974842767*lx2 - 22.931769722814497*lx3)
           + x3*(-156.05371972800287 + 15.715759007631444*lx - 52.68176100628931*lx2 - 9.239558707643814*lx3)
           + x2*(657.2565997888067 + 236.1380281690141*lx + 20.77998017839445*lx2 - 1.253121452894438*lx3)
           + x*(187.80684212073498 - 153.01811971500075*lx - 1.3333333333333333*lx2 - 0.2962962962962963*lx3)
           + x*(53.07484324195591 - 20.*lx - 6.444444444444445*lx2 + 6.222222222222222*lx3)
           + x3*(-386.1741472172352 + 663.353947368421*lx - 13.659211927582534*lx2 + 47.87750556792873*lx3)
           + _nf*(-1.9753086419753085 - 1.6790123456790123*xb - 1.1642728904847397*x4 + 0.10468594217347957*x5
                  - 0.009333333333333334*x6 - 8.*lx + 0.4444444444444444*lx2
                  + x*(-3.555559830375158 - 0.8888888888888888*lx + 0.4444444444444444*lx2)
                  + x2*(-0.4012345679012346 - 8.9925705794948*lx + 0.8817787418655098*lx2)
                  + x3*(10.359047619047619 - 4.647992530345472*lx + 1.703971119133574*lx2));
  }
  double C2nspff::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C2nspff::Local(double const& x) const
  {
    const double A1 = -30.608535416484305 + 3.4489756184793317*_nf;
    return A1 + _A2 * log(1 - x);
  }

  //_________________________________________________________________________________
  C2nsmff::C2nsmff(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = 14.926669450170849 + 5.530864197530864*_nf;
  }
  double C2nsmff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double x4   = x * x3;
    const double x5   = x * x4;
    const double x6   = x * x5;
    const double xb   = 1 - x;
    const double xb2  = xb * xb;
    const double xb3  = xb * xb2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    return -256.09239999688725 + 375.93360231139445*x4 - 29.662873954228793*x5 + 2.4864127231800848*x6
           - 22.22222222222222*l1x + 7.111111111111111*l1x2
           + xb*(141.04140553909332 + 47.111111111111114*l1x - 3.5555555555555554*l1x2)
           + xb3*(3.9346330275229358*l1x - 0.7608247422680412*l1x2)
           + xb2*(-12.865738999623918*l1x + 2.426546391752577*l1x2)
           + 72.31521361833259*lx - 10.88888888888889*lx2 - 3.8518518518518516*lx3
           + x2*(56.065623348965666 - 30.13516702067784*lx - 7.172955974842767*lx2 - 22.931769722814497*lx3)
           + x3*(-156.05371972800287 + 15.715759007631444*lx - 52.68176100628931*lx2 - 9.239558707643814*lx3)
           - x2*(657.2565997888067 + 236.1380281690141*lx + 20.77998017839445*lx2 - 1.253121452894438*lx3)
           + x*(187.80684212073498 - 153.01811971500075*lx - 1.3333333333333333*lx2 - 0.2962962962962963*lx3)
           - x*(53.07484324195591 - 20.*lx - 6.444444444444445*lx2 + 6.222222222222222*lx3)
           - x3*(-386.1741472172352 + 663.353947368421*lx - 13.659211927582534*lx2 + 47.87750556792873*lx3)
           + _nf*(-1.9753086419753085 - 1.6790123456790123*xb - 1.1642728904847397*x4 + 0.10468594217347957*x5
                  - 0.009333333333333334*x6 - 8.*lx + 0.4444444444444444*lx2
                  + x*(-3.555559830375158 - 0.8888888888888888*lx + 0.4444444444444444*lx2)
                  + x2*(-0.4012345679012346 - 8.9925705794948*lx + 0.8817787418655098*lx2)
                  + x3*(10.359047619047619 - 4.647992530345472*lx + 1.703971119133574*lx2));
  }
  double C2nsmff::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C2nsmff::Local(double const& x) const
  {
    const double A1 = -30.608535416484305 + 3.4489756184793317*_nf;
    return A1 + _A2 * log(1 - x);
  }

  //_________________________________________________________________________________
  C2psff::C2psff():
    Expression()
  {
  }
  double C2psff::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double x3  = x * x2;
    const double x4  = x * x3;
    const double x5  = x * x4;
    const double x6  = x * x5;
    const double lx  = log(x);
    const double lx2 = lx * lx;
    const double lx3 = lx * lx2;
    return -94.87929671304775 - 0.17378497790868924*x4 + 0.0074487895716946*x5 - 0.000427715996578272*x6
           - 96.*lx - 14.666666666666666*lx2 + 9.777777777777779*lx3
           + (2*(-1.4599875154038369 + 3.5555555555555554*lx + 7.111111111111111*lx2))/x
           + 2*x2*(11.416978776529339 - 5.589016829052259*lx - 5.3288241415192505*lx2 + 0.000177367860943597*lx3)
           + 2*x3*(3.4597249508840866 - 1.0824053452115814*lx + 0.12507896399241947*lx2 + 0.024505183788878417*lx3)
           + 2*x*(34.10631502319054 - 18.666666666666668*lx - 7.333333333333333*lx2 + 4.888888888888889*lx3);
  }

  //_________________________________________________________________________________
  C2qgff::C2qgff():
    Expression()
  {
  }
  double C2qgff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double x4   = x * x3;
    const double x5   = x * x4;
    const double x6   = x * x5;
    const double xb   = 1 - x;
    const double xb2  = xb * xb;
    const double xb3  = xb * xb2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x * l1x2;
    return 2*(1290.4934065822974 + 53.05125148986889*x4 - 3.6024653312788906*x5 - 0.587737843551797*x6
              - 48.30935289373048*l1x - 4.444444444444445*l1x2 + 1.4814814814814814*l1x3
              + xb3*(-86.93528505392912*l1x + 37.6297520661157*l1x2 - 3.840566037735849*l1x3)
              + xb*(-441.62535053517735 - 33.19824178261937*l1x + 21.333333333333332*l1x2 + 1.4814814814814814*l1x3)
              + xb2*(-66.36817653890824*l1x + 23.107591754650578*l1x2 + 2.869815668202765*l1x3)
              + 381.3333333333333*lx + 21.333333333333332*lx2 - 45.629629629629626*lx3
              + (-152.7222196816283 + 133.33333333333337*lx - 282.6666666666667*lx2 - 106.66666666666667*lx3)/x
              + x*(-1490.8212970756279 + 139.55555555555554*lx + 84.88888888888889*lx2 - 89.18518518518519*lx3)
              + x3*(15.196861626248216 - 110.80171184022825*lx - 22.529929577464788*lx2 - 5.923846153846154*lx3)
              + x2*(304.1022632020117 - 67.32855680655067*lx - 5.704938271604938*lx2 - 0.2122492080253432*lx3));
  }

  //_________________________________________________________________________________
  C2gqff::C2gqff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C2gqff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double x4   = x * x3;
    const double x5   = x * x4;
    const double x6   = x * x5;
    const double xb   = 1 - x;
    const double xb2  = xb * xb;
    const double xb3  = xb * xb2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    const double l1x3 = l1x * l1x2;
    return -33.57010886514981 + 33.326492537313435*x4 - 0.5029126213592233*x5 - 0.7086999022482894*x6
           + 1.7826740018155984*l1x - 3.5*l1x2 - 0.5555555555555556*l1x3
           + xb2*(41.88318144159072*l1x - 5.209855564995752*l1x2 - 1.0109140518417463*l1x3)
           + xb*(-79.15825323175919 - 17.232014670297865*l1x + 8.333333333333334*l1x2 + 1.1111111111111112*l1x3)
           + xb3*(13.599318955732123*l1x - 3.710654936461388*l1x2 + 3.189528795811518*l1x3)
           - 107.16446889494088*lx - 24.333333333333332*lx2 + 8.555555555555555*lx3
           + (-3.2849719096586334 + 8.*lx + 16.*lx2)/x
           + x2*(-179.78178152641595 - 52.541380092890414*lx - 45.45358401880141*lx2 - 9.994087213599409*lx3)
           + x3*(226.18277680140596 - 148.6697435897436*lx - 12.601598173515981*lx2 - 6.588960342979636*lx3)
           + x*(-20.762874884364063 - 127.67106221011824*lx + 2.6666666666666665*lx2 + 66.88888888888889*lx3)
           + _nf*(-6.0918268703926 - 0.2945765230312036*x4 + 0.04507042253521127*x5 - 0.015025041736227046*x6
                  + 1.1111111111111112*l1x + 0.3333333333333333*l1x2
                  + xb*(4.876032563689201 - 1.5555555555555556*l1x - 0.6666666666666666*l1x2)
                  + xb3*(2.042*l1x - 0.17929759704251386*l1x2) + xb2*(0.6025934401220442*l1x + 0.6704805491990846*l1x2)
                  - 3.3333333333333335*lx + 0.3333333333333333*lx2
                  + x*(7.5077529937967284 + 7.333333333333333*lx - 0.6666666666666666*lx2)
                  + x3*(0.5140018066847335 + 2.0420121614151463*lx - 0.17932148626817448*lx2)
                  + x2*(-2.8811926605504588 - 8.286269430051814*lx + 0.6704834605597965*lx2));
  }

  //_________________________________________________________________________________
  C2ggff::C2ggff(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = 33.585006262884406 + 12.444444444444445*_nf;
  }
  double C2ggff::Regular(double const& x) const
  {
    const double x2   = x * x;
    const double x3   = x * x2;
    const double x4   = x * x3;
    const double x5   = x * x4;
    const double x6   = x * x5;
    const double xb   = 1 - x;
    const double xb2  = xb * xb;
    const double xb3  = xb * xb2;
    const double lx   = log(x);
    const double lx2  = lx * lx;
    const double lx3  = lx * lx2;
    const double l1x  = log(1 - x);
    const double l1x2 = l1x * l1x;
    return 1795.4864691973607 + 5767.008457081538*x4 - 417.0504145663607*x5 + 25.780639491607186*x6 - 6.*l1x + 36.*l1x2
           + xb*(444.8228917453772 + 216.*l1x - 18.*l1x2) + xb3*(287.32341526520054*l1x + 0.022598870056497175*l1x2)
           + xb2*(5.533766233766234*l1x + 66.38206627680312*l1x2) + 1411.1762640653615*lx + 33.*lx2 - 132.*lx3
           + x3*(4684.30480352509 - 11480.558892997977*lx - 103.32680320569902*lx2 - 967.0959537572254*lx3)
           + x*(-2387.244708485864 + 120.64747186927704*lx - 75.*lx2 - 744.*lx3)
           + (-296.95832761699717 - 56.17626406536148*lx - 660.*lx2 - 240.*lx3)/x
           + x2*(-9721.782429932475 - 3592.3978209442635*lx - 591.4944707740916*lx2 - 21.475138121546962*lx3)
           + _nf*(-127.7777777777778 + 37.55555555555556*xb - 3.928380187416332*x4 + 0.15381831342743396*x5
                  - 0.007829977628635347*x6 + 2.*l1x - 0.000020835503698301906*xb2*l1x - 0.00022306491188935982*xb3*l1x
                  - 50.*lx + 7.333333333333333*lx2 + 9.777777777777779*lx3
                  + (3.432098765432099 - 27.11111111111111*lx - 0.8888888888888888*lx2)/x
                  + x2*(49.03990326481257 - 25.708160442600278*lx + 7.786840301576422*lx2 + 0.029953917050691243*lx3)
                  + x3*(0.643843498273878 + 15.898224852071007*lx - 4.557580174927113*lx2 + 1.7672833495618305*lx3)
                  + x*(83.33321256038647 - 78.*lx - 19.333333333333332*lx2 + 9.777777777777779*lx3));
  }
  double C2ggff::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C2ggff::Local(double const& x) const
  {
    const double A1 = -62.104684476395086 + 7.760195141578497*_nf;
    return A1 + _A2 * log(1 - x);
  }

  //_________________________________________________________________________________
  C3nspff::C3nspff(int const& nf):
    Expression(),
    _nf(nf)
  {
    _A2 = 0;
  }
  double C3nspff::Regular(double const&) const
  {
    return 0;
  }
  double C3nspff::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C3nspff::Local(double const& x) const
  {
    const double ln1mx = log( 1 - x );
    const double A1    = 0;

    return A1 + _A2 * ln1mx;
  }

  //_________________________________________________________________________________
  C3nsmff::C3nsmff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3nsmff::Regular(double const&) const
  {
    return 0;
  }
  double C3nsmff::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C3nsmff::Local(double const& x) const
  {
    const double ln1mx = log( 1 - x );
    const double A1    = 0;

    return A1 + _A2 * ln1mx;
  }

  //_________________________________________________________________________________
  C3psff::C3psff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3psff::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C3pvff::C3pvff():
    Expression()
  {
  }
  double C3pvff::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C3qgff::C3qgff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3qgff::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C3gqff::C3gqff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3gqff::Regular(double const&) const
  {
    return 0;
  }

  //_________________________________________________________________________________
  C3ggff::C3ggff(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double C3ggff::Regular(double const&) const
  {
    return 0;
  }
  double C3ggff::Singular(double const& x) const
  {
    return _A2 / ( 1 - x );
  }
  double C3ggff::Local(double const& x) const
  {
    const double ln1mx = log( 1 - x );
    const double A1    = 0;

    return A1 + _A2 * ln1mx;
  }
}
