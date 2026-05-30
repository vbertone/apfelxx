//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionsunp_sl_krk.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  // ============================================================================
  // P1qiqiKrk
  // ============================================================================
  // Source: Pqiqi[Krk][2] from PrecomputedSplittingFunctions.m
  //
  // Distributions (before x4 factor):
  //   D0 coeff: CF*nf*TR*(-20/9) + CA*CF*(67/9 - pi^2/3)
  //   D1 coeff: CF*nf*TR*(8/3)   + CA*CF*(-22/3)
  //   delta:    CF*nf*TR*(-3-4*pi^2/9) + CA*CF*(17/2+11*pi^2/9-3*zeta3)
  //             + CF^2*(3/8-pi^2/2+6*zeta3)
  //
  // Singular(x) = (_a0 + _a1*log(1-x)) / (1-x)
  // Local(x)    = _delta + _a0*log(1-x) + (_a1/2)*log^2(1-x)
  //_________________________________________________________________________________
  P1qiqiKrk::P1qiqiKrk(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a0    = 4. * ( CF * _nf * TR * ( - 20. / 9. )
                    + CA * CF * ( 67. / 9. - Pi2 / 3. ) );
    _a1    = 4. * ( CF * _nf * TR * ( 8. / 3. )
                    + CA * CF * ( - 22. / 3. ) );
    _delta = 4. * ( CF * _nf * TR * ( - 3. - 4. * Pi2 / 9. )
                    + CA * CF * ( 17. / 2. + 11. * Pi2 / 9. - 3. * zeta3 )
                    + CF * CF * ( 3. / 8. - Pi2 / 2. + 6. * zeta3 ) );
  }
  double P1qiqiKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    return 4. * (
             // CF*TR term (= Pqiqk regular part)
             CF * TR * ( - 5. + 26. / ( 9. * x ) + 5. * x - 26. * x2 / 9.
                         + ( 4. / ( 3. * x ) + 2. * x + 4. * x2 / 3. ) * lnx
                         + ( - 1. - x ) * lnx2 )
             // CF*nf*TR regular terms
             + CF * _nf * TR * ( 4. * ( 1. + 4. * x ) / 9.
                                 - 4. * ( 1. + x ) * ln1mx / 3.
                                 - 4. * ( 1. + x2 ) * lnx / ( 3. * ( 1. - x ) ) )
             // CA*CF regular terms
             + CA * CF * ( ( 10. - 77. * x ) / 9. + Pi2 * ( 1. + x ) / 6.
                           + 11. * ( 1. + x ) * ln1mx / 3.
                           + ( 14. + 8. * x2 ) * lnx / ( 3. * ( 1. - x ) )
                           + ( 1. + x2 ) * lnx2 / ( 2. * ( 1. - x ) ) )
             // CF^2 regular terms
             + CF * CF * ( 5. * ( x - 1. )
                           - ( 3. + 2. * x - 2. * x2 ) * lnx / ( 1. - x )
                           - 2. * ( 1. + x2 ) * ln1mx * lnx / ( 1. - x )
                           + ( - 1. - x ) * lnx2 / 2. )
           );
  }
  double P1qiqiKrk::Singular(double const& x) const
  {
    return ( _a0 + _a1 * log(1 - x) ) / ( 1 - x );
  }
  double P1qiqiKrk::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return _delta + _a0 * ln1mx + ( _a1 / 2. ) * ln1mx * ln1mx;
  }

  // ============================================================================
  // P1qiqkKrk
  // ============================================================================
  // Source: Pqiqk[Krk][2] = Pqiqbk[Krk][2]
  // Purely regular, no nf dependence.
  //_________________________________________________________________________________
  P1qiqkKrk::P1qiqkKrk():
    Expression()
  {
  }
  double P1qiqkKrk::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double lnx = log(x);
    return 4. * CF * TR * ( - 5. + 26. / ( 9. * x ) + 5. * x - 26. * x2 / 9.
                            + ( 4. / ( 3. * x ) + 2. * x + 4. * x2 / 3. ) * lnx
                            + ( - 1. - x ) * lnx * lnx );
  }

  // ============================================================================
  // P1qiqbiKrk
  // ============================================================================
  // Source: Pqiqbi[Krk][2] = Pqiqk part + CF^2 part + CA*CF part.
  // Purely regular, no nf dependence.
  //
  // The CF^2 and CA*CF contributions combine to:
  //   (2*CF^2 - CA*CF) * [2*(1-x) + (1+x)*lnx + pqqmxc*S2x]
  // where pqqmxc = (1+x^2)/(1+x) and S2x = standard Catani-Seymour function.
  //_________________________________________________________________________________
  P1qiqbiKrk::P1qiqbiKrk():
    Expression()
  {
  }
  double P1qiqbiKrk::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1px  = log(1 + x);
    const double pqqmxc = ( 1. + x2 ) / ( 1. + x );
    // S2x = -2*Li2(-x) + lnx^2/2 - 2*lnx*log(1+x) - pi^2/6
    const double S2x    = - 2. * dilog(-x) + lnx2 / 2. - 2. * lnx * ln1px - Pi2 / 6.;
    return 4. * (
             // CF*TR part (identical to P1qiqkKrk)
             CF * TR * ( - 5. + 26. / ( 9. * x ) + 5. * x - 26. * x2 / 9.
                         + ( 4. / ( 3. * x ) + 2. * x + 4. * x2 / 3. ) * lnx
                         + ( - 1. - x ) * lnx2 )
             // CF^2 part
             + CF * CF * ( 4. * ( 1. - x ) + 2. * ( 1. + x ) * lnx + 2. * pqqmxc * S2x )
             // CA*CF part
             + CA * CF * ( - 2. * ( 1. - x ) - ( 1. + x ) * lnx - pqqmxc * S2x )
           );
  }

  // ============================================================================
  // P1qigKrk
  // ============================================================================
  // Source: Pqig[Krk][2]. Purely regular.
  // Defines pqg = (1-x)^2 + x^2, pqgmx = x^2 + (1+x)^2.
  //_________________________________________________________________________________
  P1qigKrk::P1qigKrk(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1qigKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double ln1px = log(1 + x);
    const double pqg   = ( 1. - x ) * ( 1. - x ) + x2;
    const double pqgmx = x2 + ( 1. + x ) * ( 1. + x );
    return 4. * (
             // nf*TR^2 term
             - 59. * _nf * TR * TR * pqg / 36.
             // CF*TR terms
             + CF * TR * ( 7. / 4. - 8. * x + 9. * x2 / 2.
                           - 3. * pqg * ln1mx
                           - pqg * ln1mx * ln1mx
                           + ( 2. - 8. * x + 6. * x2 ) * lnx
                           + ( 1. / 2. - x + 2. * x2 ) * lnx2 )
             // CA*TR terms
             + CA * TR * ( Pi2 * ( - 1. - 2. * x2 ) / 3.
                           + ( 208. - 91. * x + 686. * x2 - 390. * x3 ) / ( 72. * x )
                           + pqg * ln1mx * ln1mx
                           + ( 4. / ( 3. * x ) + 4. * x + 31. * x2 / 3. ) * lnx
                           - 2. * pqg * ln1mx * lnx
                           + ( - 1. - 2. * x ) * lnx2
                           - 2. * pqgmx * lnx * ln1px
                           - 2. * pqgmx * dilog(-x) )
           );
  }

  // ============================================================================
  // P1gqKrk
  // ============================================================================
  // Source: Pgq[Krk][2]. Purely regular.
  // Defines pgq = (1+(1-x)^2)/x, pgqmx = (1+(1+x)^2)/x.
  //_________________________________________________________________________________
  P1gqKrk::P1gqKrk(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1gqKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double ln1px = log(1 + x);
    const double pgq   = ( 1. + ( 1. - x ) * ( 1. - x ) ) / x;
    const double pgqmx = ( 1. + ( 1. + x ) * ( 1. + x ) ) / x;
    return 4. * (
             // CF*nf*TR terms
             CF * _nf * TR * ( - 7. * pgq / 12.
                               + 4. * pgq * ln1mx / 3.
                               - 4. * pgq * lnx / 3. )
             // CF^2 terms
             + CF * CF * ( - 7. + 11. / ( 2. * x ) - Pi2 * pgq / 3. + 5. * x / 4.
                           + pgq * ln1mx * ln1mx
                           - 2. * pgq * ln1mx * lnx
                           + ( - 2. + x ) * lnx2 / 2. )
             // CA*CF terms
             + CA * CF * ( Pi2 * ( 2. + x2 ) / ( 3. * x )
                           + ( - 538. + 690. * x - 93. * x2 + 208. * x3 ) / ( 72. * x )
                           - 11. * pgq * ln1mx / 3.
                           - pgq * ln1mx * ln1mx
                           + ( - 46. + 9. / x + 5. * x - 4. * x2 ) * lnx / 3.
                           + ( 2. + x ) * lnx2
                           + 2. * pgqmx * lnx * ln1px
                           + 2. * pgqmx * dilog(-x) )
           );
  }

  // ============================================================================
  // P1ggKrk
  // ============================================================================
  // Source: Pgg[Krk][2].
  //
  // Distributions (before x4 factor):
  //   D0 coeff: CA*nf*TR*(-20/9)  + CA^2*(67/9 - pi^2/3)
  //   D1 coeff: CA*nf*TR*(8/3)    + CA^2*(-22/3)
  //   delta:    CF*nf*TR*(-1) + nf^2*TR^2*(59/54)
  //             + CA*nf*TR*(-1619/216 - 2*pi^2/9)
  //             + CA^2*(4903/432 + 11*pi^2/18 + 3*zeta3)
  //_________________________________________________________________________________
  P1ggKrk::P1ggKrk(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a0    = 4. * ( CA * _nf * TR * ( - 20. / 9. )
                    + CA * CA * ( 67. / 9. - Pi2 / 3. ) );
    _a1    = 4. * ( CA * _nf * TR * ( 8. / 3. )
                    + CA * CA * ( - 22. / 3. ) );
    _delta = 4. * ( CF * _nf * TR * ( - 1. )
                    + _nf * _nf * TR * TR * ( 59. / 54. )
                    + CA * _nf * TR * ( - 1619. / 216. - 2. * Pi2 / 9. )
                    + CA * CA * ( 4903. / 432. + 11. * Pi2 / 18. + 3. * zeta3 ) );
  }
  double P1ggKrk::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double x4    = x * x3;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double ln1px = log(1 + x);
    // pgg includes the 1/(1-x) piece; terms multiplied by lnx or lnx*ln1mx
    // are all integrable at x=1 (lnx -> 0 cancels the pole).
    const double pgg   = 1. / ( 1. - x ) + 1. / x - 2. + x * ( 1. - x );
    const double pggmx = 1. / ( 1. + x ) - 1. / x - 2. - x * ( 1. + x );
    return 4. * (
             // CF*nf*TR terms
             CF * _nf * TR * ( 10. * ( x - 1. )
                               + ( - 4. - 8. / ( 3. * x ) - 4. * x + 8. * x2 / 3. ) * lnx
                               - 2. * ( 1. + x ) * lnx2 )
             // CA*nf*TR terms (regular, distributions excluded)
             + CA * _nf * TR * ( 2. * ( - 23. + 29. * x - 19. * x2 + 23. * x3 ) / ( 9. * x )
                                 - 8. * ( x3 - x2 + 2. * x - 1. ) * ln1mx / ( 3. * x )
                                 - 4. * ( 1. - x + 3. * x2 - 3. * x3 + x4 ) * lnx
                                 / ( 3. * ( 1. - x ) * x ) )
             // CA^2 terms (regular, distributions excluded)
             + CA * CA * ( ( - 25. - 109. * x ) / 18.
                           + Pi2 * ( 3. + 4. * x + 2. * x2 + 2. * x3 ) / ( 3. * ( 1. + x ) )
                           + 22. * ( x3 - x2 + 2. * x - 1. ) * ln1mx / ( 3. * x )
                           + ( 11. - 47. * x + 69. * x2 - 77. * x3 + 55. * x4 ) * lnx
                           / ( 3. * x * ( 1. - x ) )
                           - 4. * pgg * ln1mx * lnx
                           - 2. * ( 1. + x - x2 ) * ( 1. + x - x2 ) * lnx2 / ( x2 - 1. )
                           - 4. * pggmx * lnx * ln1px
                           - 4. * pggmx * dilog(-x) )
           );
  }
  double P1ggKrk::Singular(double const& x) const
  {
    return ( _a0 + _a1 * log(1 - x) ) / ( 1 - x );
  }
  double P1ggKrk::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return _delta + _a0 * ln1mx + ( _a1 / 2. ) * ln1mx * ln1mx;
  }
}
