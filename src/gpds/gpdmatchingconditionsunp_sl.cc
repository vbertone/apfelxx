//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdmatchingconditionsunp_sl.h"
#include "apfel/constants.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  AS1gpdgH_0::AS1gpdgH_0(double const& xi):
    Expression(1 / xi)
  {
  }
  double AS1gpdgH_0::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;

    // Parton-level skewness. Only computed for s < 1.
    const double s = kappa * y;
    if (s >= 1 || (kappa < 1 && y >= 1))
      return 0;

    // Small kappa: series in s, which avoids the 1/s cancellations
    if (kappa < 1e-3)
      {
        const double y2 = y * y;
        const double y3 = y2 * y;
        const double y4 = y2 * y2;
        const double y5 = y4 * y;
        const double s2 = s * s;
        const double ly = log(y);
        return CF * ( 2 * ( 1 + ( 1 - y ) * ( 1 - y ) ) * ( - 1 - 2 * ly ) / y
                      - 4 * ( y - 1 ) * ( - 1 + 8 * y - 10 * y2 + 6 * y3 + 3 * ( y - 1 ) * y2 * ly ) / 3 / y3 * s2
                      - ( y - 1 ) * ( - 2 - 32 * y + 133 * y2 - 7 * y3 - 187 * y4 + 155 * y5 + 60 * ( y - 1 ) * y4 * ly ) / 15 / y5 * s2 * s2 );
      }

    // Removable singularity at y = s^2 (kappa > 1): average over two
    // nearby points
    if (kappa > 1 && std::abs(1 - kappa * s) < 1e-5)
      return ( Regular(y * ( 1 + 2e-5 )) + Regular(y * ( 1 - 2e-5 )) ) / 2;

    // The HDR result in parton-in-parton variables (x, t) is f(x, xi)
    // Theta(xi) + f(x, -xi) Theta(-xi). For kappa < 1 (DGLAP) this
    // gives f(y, s) + f(y, -s), for kappa > 1 (ERBL) the combination
    // M(y) - M(-y) gives f(y, s) - f(-y, s).
    double res = 0;
    for (int i = 0; i < 2; i++)
      {
        const double x    = (i == 1 && kappa > 1 ? - y : y);
        const double t    = (i == 1 && kappa < 1 ? - s : s);
        const double t2   = t * t;
        const double d    = x - t2;
        const double x2t2 = x * x - t2;
        res += (i == 1 && kappa > 1 ? - 1 : 1)
               * ( x2t2 - 2 * ( x + t ) * ( x + t ) * ( 1 - t2 - t * ( 1 - x ) ) / d * log(( x + t ) / ( 1 + t ))
                   - ( 1 + t ) * ( 1 - x ) * x2t2 / d * log(std::abs(x2t2 / ( 1 - t2 ))) ) / x / t / ( 1 + t );
      }
    return CF * res;
  }

  //_________________________________________________________________________________
  AS1gpdHH_0::AS1gpdHH_0(double const& xi):
    Expression(1 / xi)
  {
  }
  double AS1gpdHH_0::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;

    // Parton-level skewness. Only computed for s < 1.
    const double s = kappa * y;
    if (s >= 1 || (kappa < 1 && y >= 1))
      return 0;

    // Small s (DGLAP): series in s, which avoids the 1/s cancellations
    if (kappa < 1 && s < 1e-3)
      {
        const double y2 = y * y;
        const double s2 = s * s;
        const double l1 = log(1 - y);
        return CF * ( 2 * ( 1 + y2 ) * ( - 1 - 2 * l1 ) / ( 1 - y )
                      + 4 * ( 3 - 2 * y + 6 * y2 + 3 * ( y2 - 1 ) * l1 ) / 3 / ( y - 1 ) * s2
                      + ( - 45 - 48 * y + 155 * y2 + 60 * ( y2 - 1 ) * l1 ) / 15 / ( y - 1 ) * s2 * s2 ) - Singular(y);
      }

    // The HDR result in parton-in-parton variables (x, t) is f(x, xi)
    // Theta(xi) + f(x, -xi) Theta(-xi), plus-prescripted. For kappa <
    // 1 (DGLAP) this gives f(y, s) + f(y, -s), for kappa > 1 (ERBL)
    // the combination M(y) - M(-y) gives f(y, s) - f(-y, s).
    double res = 0;
    for (int i = 0; i < 2; i++)
      {
        const double x   = (i == 1 && kappa > 1 ? - y : y);
        const double t   = (i == 1 && kappa < 1 ? - s : s);
        const double t2  = t * t;
        const double omx = 1 - x;
        res += (i == 1 && kappa > 1 ? - 1 : 1)
               * ( 2 * ( x - t ) * ( 1 + t - t * ( x + t ) ) / t / omx / ( 1 + t ) * log(omx / ( 1 + t ))
                   - ( 1 + t2 - omx * omx ) / t / omx * log(( 1 - x * x ) / ( 1 - t2 ))
                   - ( x + t ) / t / ( 1 + t ) );
      }
    return CF * res - Singular(y);
  }
  double AS1gpdHH_0::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa >= 1 || y >= 1)
      return 0;

    // Singular part at y -> 1: (- 8 ln(1 - y) + B) / (1 - y)
    const double B = (kappa < 1e-8 ? - 4 : 2 * ( ( 1 + kappa ) * ( 1 + kappa ) * log1p(- kappa) - ( 1 - kappa ) * ( 1 - kappa ) * log1p(kappa) ) / kappa);
    return CF * ( - 8 * log(1 - y) + B ) / ( 1 - y );
  }
  double AS1gpdHH_0::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa >= 1)
      return 0;

    // Coefficient B of the singular part and coefficient E of
    // delta(1 - y). E combines the HDR delta(1 - x) term with the
    // conversion of the plus prescription at fixed parton-level
    // skewness to the one at fixed kappa.
    const double lm = log1p(- kappa);
    const double lp = log1p(kappa);
    const double B  = (kappa < 1e-8 ? - 4 : 2 * ( ( 1 + kappa ) * ( 1 + kappa ) * lm - ( 1 - kappa ) * ( 1 - kappa ) * lp ) / kappa);
    const double E  = (kappa < 1e-8 ? 4 : 4 - ( ( 1 + kappa ) * ( 1 + kappa ) * lm * lm - ( 1 - kappa ) * ( 1 - kappa ) * lp * lp ) / kappa);
    const double l1 = log(1 - y);
    return CF * ( E - 4 * l1 * l1 + B * l1 );
  }

  //_________________________________________________________________________________
  AS1gpdHH_L::AS1gpdHH_L(double const& xi):
    Pgpd0qq(xi)
  {
  }
  double AS1gpdHH_L::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa >= 1)
      return Pgpd0qq::Local(y);

    // Cusp term 2 C_F (chi coth(chi) - 1) delta(1 - y), with chi =
    // ln((1 + kappa) / (1 - kappa))
    const double chi = log(( 1 + kappa ) / ( 1 - kappa ));
    return Pgpd0qq::Local(y) + 2 * CF * (chi < 1e-4 ? chi * chi / 3 : chi / tanh(chi) - 1);
  }
}
