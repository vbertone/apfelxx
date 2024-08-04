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
      return 2 * TR / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1qg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR / pow(kappa, 2) * log(1 - kappa * y);
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
      return (y <= 1 ? - 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - CA * ( 1 + kappa * ( y + kappa * ( - 3 + ky ) ) ) / kappa / pow(1 + ky, 2) : - 2 * CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 3 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2));
  }
  double Cgtmd1gg::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1gg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1gg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1polns::Cgtmd1polns(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polns::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1polns::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1polqq::Cgtmd1polqq(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polqq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1polqq::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1polqg::Cgtmd1polqg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polqg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * TR * ( 3 + ky ) / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1polqg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1polqg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1polgq::Cgtmd1polgq(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polgq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 4 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1polgg::Cgtmd1polgg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polgg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * CA * ( 3 + ky ) / pow(1 + ky, 2) : - 8 * CA * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1polgg::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1polgg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1polgg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1lingg::Cgtmd1lingg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1lingg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * CA * kappa * ( 1 - ky ) / pow(1 + ky, 2) : - 8 * CA * ky2 * ( 1 - kappa ) / pow(1 - ky2, 2));
  }
  double Cgtmd1lingg::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1lingg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA * kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1lingg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1linunpgq::Cgtmd1linunpgq(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linunpgq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF * ( 1 - y ) / y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 4 * CF * ( 1 - kappa ) / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1linunpgg::Cgtmd1linunpgg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linunpgg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 4 * CA * ( 1 - y ) * ( 1 + ky2 ) / y / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * CA * ( 2 + ky * ( 1 + ky ) ) / y / pow(1 + ky, 2) : - 4 * CA * ( 1 - kappa ) * ( 1 + ky2 ) / pow(1 - ky2, 2));
  }
  double Cgtmd1linunpgg::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1linunpgg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA * kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1linunpgg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1linpolgq::Cgtmd1linpolgq(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linpolgq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 - y ) * kappa / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * ( 1 - kappa ) / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1linpolgg::Cgtmd1linpolgg(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linpolgg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 4 * CA * kappa * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - CA * ( 3 + ky ) / pow(1 + ky, 2) : - 4 * CA * ( 1 - kappa ) / pow(1 - ky2, 2));
  }
  double Cgtmd1linpolgg::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1linpolgg::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - CA / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1linpolgg::LocalPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - CA * log(1 - kappa * y) / kappa;
    else
      return 0;
  }
}
