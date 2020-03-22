//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdevolutionkernels.h"
#include "apfel/constants.h"

#include <iostream>

namespace apfel
{
  //_________________________________________________________________________________
  Pgpd0nsDGLAP::Pgpd0nsDGLAP(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Pgpd0nsDGLAP::Singular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      return 2 * CF * ( 2 / ( 1 - x )
                        + ( 1 - kappa ) / ( 2 * kappa * ( 1 + kappa * x ) )
                        - ( 1 + kappa ) / ( 2 * kappa * ( 1 - kappa * x ) ) );
  }
  double Pgpd0nsDGLAP::Local(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      return 4 * CF * ( log( 1 - x ) + ( - ( 1 + kappa ) * log( 1 - kappa * x )
                                         - ( 1 - kappa ) * log( 1 + kappa * x )
                                       ) / pow(2 * kappa, 2) );
  }

  //_________________________________________________________________________________
  Pgpd0nsERBL::Pgpd0nsERBL(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Pgpd0nsERBL::Singular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 2 * CF * ( 2 / ( 1 - x )
                        + ( 1 - kappa ) / ( 2 * kappa * ( 1 + kappa * x ) ) );
    else
      return 0;
  }

  //_________________________________________________________________________________
  Pgpd0qgDGLAP::Pgpd0qgDGLAP(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Pgpd0qgDGLAP::Regular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      {
        const double fact  = 1 - pow(kappa * x, 2);
        return 4 * TR * ( 1 - 2 * x * ( 1 - x ) / fact ) / fact;
      }
  }

  //_________________________________________________________________________________
  Pgpd0qgERBL::Pgpd0qgERBL(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Pgpd0qgERBL::Regular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 2 * TR * ( 1 + kappa ) * ( - 1 + kappa + ( - 2 + kappa ) * kappa * x ) /( pow(kappa, 3) * x * pow(1 + kappa * x, 2) );
    else
      return 0;
  }

  //_________________________________________________________________________________
  Pgpd0gqDGLAP::Pgpd0gqDGLAP(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Pgpd0gqDGLAP::Regular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      return 2 * CF * ( pow(1 - x, 2) / ( 1 - pow(kappa * x, 2) ) + 1 ) /x;
  }

  //_________________________________________________________________________________
  Pgpd0gqERBL::Pgpd0gqERBL(double const& xi):
    Expression(),
    _xi(xi)
  {
  }
  double Pgpd0gqERBL::Regular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return CF * ( 2 * kappa + ( - 1 + ( - 2 + kappa ) * kappa ) * x ) / ( kappa * x * ( 1 + kappa * x ) );
    else
      return 0;
  }

  //_________________________________________________________________________________
  Pgpd0ggDGLAP::Pgpd0ggDGLAP(int const& nf, double const& xi):
    Expression(),
    _nf(nf),
    _xi(xi)
  {
  }
  double Pgpd0ggDGLAP::Regular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      return 4 * CA / x;
  }
  double Pgpd0ggDGLAP::Singular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      {
        const double x2     = x * x;
        const double kappa2 = kappa * kappa;
        return 4 * CA * ( ( - 2 + ( 1 + kappa2 ) * x - ( 1 - kappa2 ) * x2 ) / pow(1 - kappa2 * x2, 2) + 1 / ( 1 - x ) );
      }
  }
  double Pgpd0ggDGLAP::Local(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 0;
    else
      {
        const double x2     = x * x;
        const double kappa2 = kappa * kappa;
        const double kappa3 = kappa * kappa2;
        return 4 * CA * log( 1 - x ) + 11 / 3. * CA - 4 / 3. * _nf * TR
               + (CA*(-2*(1/kappa2 + 1)*(-1 + x)
                      +(1/kappa3 + 3*kappa*x2)*log(((-1 + kappa)*(1 + kappa*x))/((1 + kappa)*(-1 + kappa*x)))
                      +1/kappa*(3 + x2)*(-log(1 - (kappa*(-1 + x))/(-1 + kappa2*x)) +
                                         log(1 + (kappa*(-1 + x))/(-1 + kappa2*x)))))
               /(-1 + kappa2*x2);
      }
  }

  //_________________________________________________________________________________
  Pgpd0ggERBL::Pgpd0ggERBL(int const& nf, double const& xi):
    Expression(),
    _nf(nf),
    _xi(xi)
  {
  }
  double Pgpd0ggERBL::Regular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 2 * CA / x
             + CA * ( - 1 + 3 * pow(kappa, 2) - kappa * ( 2 + pow(1 - kappa, 2) * kappa ) * x )
             / ( pow(kappa, 3) * x * pow(1 + kappa * x, 2) );
    else
      return 0;
  }
  double Pgpd0ggERBL::Singular(double const& x) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 2 * CA / ( 1 - x );
    else
      return 0;
  }
  double Pgpd0ggERBL::Local(double const&) const
  {
    const double kappa = _xi / _extvar;
    if (kappa > 1)
      return 11 / 3. * CA - 4 / 3. * _nf * TR;
    else
      return 0;
  }
}
