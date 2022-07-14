//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gmatchingfunctions.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  Cgtmd1ns::Cgtmd1ns(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ns::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1ns::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qq::Cgtmd1qq(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1qq::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qg::Cgtmd1qg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * TR * ( 1 - ky ) / kappa / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * pow(y, 2) / pow(1 - ky2, 2));
  }
  double Cgtmd1qg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / ( 1 - kappa * y ) / kappa;
    else
      return 0;
  }
  double Cgtmd1qg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR * log( 1 - kappa * y ) / pow(kappa, 2);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1gq::Cgtmd1gq(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - pow(kappa, 2) ) * y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * ( 1 - pow(kappa, 2) ) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1gg::Cgtmd1gg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? CA * ( 1 + kappa * ( y + kappa * ( - 5 + 3 * ky ) ) ) / 2 / kappa / pow(1 + ky, 2) : CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 7 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2));
  }
  double Cgtmd1gg::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1gg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return CA * ( 1 + 3 * pow(kappa, 2) ) / ( 1 - kappa * y ) / 2 / kappa;
    else
      return 0;
  }
  double Cgtmd1gg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return CA * ( 1 + 3 * pow(kappa, 2) ) * log( 1 - kappa * y ) / 2 / pow(kappa, 2);
    else
      return 0;
  }
}
