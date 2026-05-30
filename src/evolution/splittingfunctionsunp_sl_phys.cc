//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionsunp_sl_phys.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  // ============================================================================
  // P1qiqiPHYS
  // ============================================================================
  // Source: Pqiqi[PHYS][2] from PrecomputedSplittingFunctions.m
  //
  // Distributions (before x4 factor):
  //   D0 coeff:  CF*nf*TR*(-20/9) + CA*CF*(67/9 - pi^2/3)
  //   D1 coeff:  CF*nf*TR*(4/3)   + CA*CF*(-11/3)
  //   delta:     CF*nf*TR*(-2-2*pi^2/9) + CA*CF*(23/4+11*pi^2/18-3*zeta3)
  //              + CF^2*(3/8-pi^2/2+6*zeta3)
  //
  // Singular(x) = (_a0 + _a1*log(1-x)) / (1-x)
  // Local(x)    = _delta + _a0*log(1-x) + (_a1/2)*log^2(1-x)
  //_________________________________________________________________________________
  P1qiqiPHYS::P1qiqiPHYS(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a0    = 4. * ( CF * _nf * TR * ( - 20. / 9. )
                    + CA * CF * ( 67. / 9. - Pi2 / 3. ) );
    _a1    = 4. * ( CF * _nf * TR * ( 4. / 3. )
                    + CA * CF * ( - 11. / 3. ) );
    _delta = 4. * ( CF * _nf * TR * ( - 2. - 2. * Pi2 / 9. )
                    + CA * CF * ( 23. / 4. + 11. * Pi2 / 18. - 3. * zeta3 )
                    + CF * CF * ( 3. / 8. - Pi2 / 2. + 6. * zeta3 ) );
  }
  double P1qiqiPHYS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    return 4. * (
             // CF*TR term (= P1qiqkPHYS regular part)
             CF * TR * ( - 5. + 13. / ( 9. * x ) + 5. * x - 13. * x2 / 9.
                         + ( - 1. + 4. * x2 / 3. ) * lnx
                         + ( - 1. - x ) * lnx2 )
             // CF*nf*TR regular terms
             + CF * _nf * TR * ( 4. * ( 1. + 4. * x ) / 9.
                                 - 2. * ( 1. + x ) * ln1mx / 3.
                                 - 2. * ( 1. + x2 ) * lnx / ( 3. * ( 1. - x ) ) )
             // CA*CF regular terms
             + CA * CF * ( ( 10. - 77. * x ) / 9. + Pi2 * ( 1. + x ) / 6.
                           + 11. * ( 1. + x ) * ln1mx / 6.
                           + ( 17. + 5. * x2 ) * lnx / ( 6. * ( 1. - x ) )
                           + ( 1. + x2 ) * lnx2 / ( 2. * ( 1. - x ) ) )
             // CF^2 regular terms
             + CF * CF * ( 5. * ( x - 1. )
                           - ( 3. + 2. * x - 2. * x2 ) * lnx / ( 1. - x )
                           - 2. * ( 1. + x2 ) * ln1mx * lnx / ( 1. - x )
                           + ( - 1. - x ) * lnx2 / 2. )
           );
  }
  double P1qiqiPHYS::Singular(double const& x) const
  {
    return ( _a0 + _a1 * log(1 - x) ) / ( 1 - x );
  }
  double P1qiqiPHYS::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return _delta + _a0 * ln1mx + ( _a1 / 2. ) * ln1mx * ln1mx;
  }

  // ============================================================================
  // P1qiqkPHYS
  // ============================================================================
  // Source: Pqiqk[PHYS][2] = Pqiqbk[PHYS][2]
  // Purely regular, no nf dependence.
  //_________________________________________________________________________________
  P1qiqkPHYS::P1qiqkPHYS():
    Expression()
  {
  }
  double P1qiqkPHYS::Regular(double const& x) const
  {
    const double x2  = x * x;
    const double lnx = log(x);
    return 4. * CF * TR * ( - 5. + 13. / ( 9. * x ) + 5. * x - 13. * x2 / 9.
                            + ( - 1. + 4. * x2 / 3. ) * lnx
                            + ( - 1. - x ) * lnx * lnx );
  }

  // ============================================================================
  // P1qiqbiPHYS
  // ============================================================================
  // Source: Pqiqbi[PHYS][2] = P1qiqkPHYS part + CF^2 part + CA*CF part.
  // The CF^2 and CA*CF parts are identical to those in P1qiqbiKrk.
  // Purely regular, no nf dependence.
  //_________________________________________________________________________________
  P1qiqbiPHYS::P1qiqbiPHYS():
    Expression()
  {
  }
  double P1qiqbiPHYS::Regular(double const& x) const
  {
    const double x2     = x * x;
    const double lnx    = log(x);
    const double lnx2   = lnx * lnx;
    const double ln1px  = log(1 + x);
    const double pqqmxc = ( 1. + x2 ) / ( 1. + x );
    // S2x = -2*Li2(-x) + lnx^2/2 - 2*lnx*log(1+x) - pi^2/6
    const double S2x    = - 2. * dilog(-x) + lnx2 / 2. - 2. * lnx * ln1px - Pi2 / 6.;
    return 4. * (
             // CF*TR part (= P1qiqkPHYS regular part)
             CF * TR * ( - 5. + 13. / ( 9. * x ) + 5. * x - 13. * x2 / 9.
                         + ( - 1. + 4. * x2 / 3. ) * lnx
                         + ( - 1. - x ) * lnx2 )
             // CF^2 part: 2*(1-x) + (1+x)*lnx + pqqmxc*S2x (times 2)
             + CF * CF * ( 4. * ( 1. - x ) + 2. * ( 1. + x ) * lnx + 2. * pqqmxc * S2x )
             // CA*CF part: -(2*(1-x) + (1+x)*lnx + pqqmxc*S2x)
             + CA * CF * ( - 2. * ( 1. - x ) - ( 1. + x ) * lnx - pqqmxc * S2x )
           );
  }

  // ============================================================================
  // P1qigPHYS
  // ============================================================================
  // Source: Pqig[PHYS][2]. Purely regular.
  // Uses dilog(x) = Li_2(x) and dilog(x*x) = Li_2(x^2).
  //_________________________________________________________________________________
  P1qigPHYS::P1qigPHYS(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1qigPHYS::Regular(double const& x) const
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
             - 29. * _nf * TR * TR * pqg / 36.
             // CF*TR terms
             + CF * TR * ( 7. / 4. - 6. * x + 4. * x2
                           - Pi2 * pqg / 3.
                           - 3. * pqg * ln1mx / 2.
                           + ( 1. / 2. - 3. * x + 3. * x2 ) * lnx
                           + ( 1. / 2. - x + 2. * x2 ) * lnx2
                           + 2. * pqg * dilog(x) )
             // CA*TR terms
             + CA * TR * ( - 2. * Pi2 * x / 3.
                           + ( 104. - 121. * x + 602. * x2 - 310. * x3 ) / ( 72. * x )
                           + ( - 1. + 31. * x2 / 3. ) * lnx
                           - 2. * pqg * ln1mx * lnx
                           + ( - 1. - 2. * x ) * lnx2
                           - 2. * pqgmx * lnx * ln1px
                           + 8. * x * dilog(x)
                           - pqgmx * dilog(x2) )
           );
  }

  // ============================================================================
  // P1gqPHYS
  // ============================================================================
  // Source: Pgq[PHYS][2]. Purely regular.
  // Uses dilog(x) = Li_2(x) and dilog(x*x) = Li_2(x^2).
  //_________________________________________________________________________________
  P1gqPHYS::P1gqPHYS(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P1gqPHYS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double ln1px = log(1 + x);
    const double pgq   = ( 1. + ( 1. - x ) * ( 1. - x ) ) / x;
    const double pgqmx = ( 1. + ( 1. + x ) * ( 1. + x ) ) / x;
    return 4. * (
             // CF*nf*TR term (no logarithms)
             CF * _nf * TR * ( - 17. * pgq / 12. )
             // CF^2 terms
             + CF * CF * ( - 3. + 1. / x + x / 4.
                           - 3. * pgq * ln1mx / 2.
                           + ( - 2. + 3. * x / 2. ) * lnx
                           - 2. * pgq * ln1mx * lnx
                           + ( - 2. + x ) * lnx2 / 2.
                           - 2. * pgq * dilog(x) )
             // CA*CF terms
             + CA * CF * ( 2. * Pi2 / 3.
                           + ( 342. - 50. / x + 9. * x + 104. * x2 ) / 72.
                           - 4. * ( 6. + x2 ) * lnx / 3.
                           + ( 2. + x ) * lnx2
                           + 2. * pgqmx * lnx * ln1px
                           - 8. * dilog(x)
                           + pgqmx * dilog(x2) )
           );
  }

  // ============================================================================
  // P1ggPHYS
  // ============================================================================
  // Source: Pgg[PHYS][2].
  //
  // Distributions (before x4 factor):
  //   D0 coeff:  CA*nf*TR*(-20/9)  + CA^2*(67/9 - pi^2/3)
  //   D1 coeff:  CA*nf*TR*(4/3)    + CA^2*(-11/3)
  //   delta:     CF*nf*TR*(-1) + nf^2*TR^2*(29/54)
  //              + CA*nf*TR*(-1013/216)
  //              + CA^2*(3385/432 + 3*zeta3)
  //_________________________________________________________________________________
  P1ggPHYS::P1ggPHYS(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a0    = 4. * ( CA * _nf * TR * ( - 20. / 9. )
                    + CA * CA * ( 67. / 9. - Pi2 / 3. ) );
    _a1    = 4. * ( CA * _nf * TR * ( 4. / 3. )
                    + CA * CA * ( - 11. / 3. ) );
    _delta = 4. * ( CF * _nf * TR * ( - 1. )
                    + _nf * _nf * TR * TR * ( 29. / 54. )
                    + CA * _nf * TR * ( - 1013. / 216. )
                    + CA * CA * ( 3385. / 432. + 3. * zeta3 ) );
  }
  double P1ggPHYS::Regular(double const& x) const
  {
    const double x2    = x * x;
    const double x3    = x * x2;
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double ln1px = log(1 + x);
    const double pgg   = 1. / ( 1. - x ) + 1. / x - 2. + x * ( 1. - x );
    const double pggmx = 1. / ( 1. + x ) - 1. / x - 2. - x * ( 1. + x );
    return 4. * (
             // CF*nf*TR terms
             CF * _nf * TR * ( - 10. + 26. / ( 9. * x ) + 10. * x - 26. * x2 / 9.
                               + ( - 2. + 8. * x2 / 3. ) * lnx
                               - 2. * ( 1. + x ) * lnx2 )
             // CA*nf*TR terms (regular, distributions excluded)
             + CA * _nf * TR * ( 2. * ( - 23. + 29. * x - 19. * x2 + 23. * x3 ) / ( 9. * x )
                                 - 4. * ( x3 - x2 + 2. * x - 1. ) * ln1mx / ( 3. * x )
                                 - 4. * ( 1. + x ) * lnx / 3. )
             // CA^2 terms (regular, distributions excluded)
             + CA * CA * ( ( - 25. - 109. * x ) / 18.
                           + Pi2 * ( 3. + 4. * x + 2. * x2 + 2. * x3 ) / ( 3. * ( 1. + x ) )
                           + 11. * ( x3 - x2 + 2. * x - 1. ) * ln1mx / ( 3. * x )
                           + ( - 25. + 11. * x - 44. * x2 ) * lnx / 3.
                           - 4. * pgg * ln1mx * lnx
                           - 2. * ( 1. + x - x2 ) * ( 1. + x - x2 ) * lnx2 / ( x2 - 1. )
                           - 4. * pggmx * lnx * ln1px
                           - 4. * pggmx * dilog(-x) )
           );
  }
  double P1ggPHYS::Singular(double const& x) const
  {
    return ( _a0 + _a1 * log(1 - x) ) / ( 1 - x );
  }
  double P1ggPHYS::Local(double const& x) const
  {
    const double ln1mx = log(1 - x);
    return _delta + _a0 * ln1mx + ( _a1 / 2. ) * ln1mx * ln1mx;
  }
}
