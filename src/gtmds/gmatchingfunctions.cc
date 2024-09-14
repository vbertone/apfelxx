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
  Cgtmd1nse::Cgtmd1nse(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1nse::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1nse::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qqe::Cgtmd1qqe(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qqe::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1qqe::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1qge::Cgtmd1qge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * TR * ( 1 - ky ) / kappa / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * pow(y, 2) / pow(1 - ky2, 2));
  }
  double Cgtmd1qge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1qge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1gqe::Cgtmd1gqe(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gqe::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - pow(kappa, 2) ) * y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * ( 1 - pow(kappa, 2) ) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1gge::Cgtmd1gge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1gge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - CA * ( 1 - 3 * pow(kappa, 2) + ( 1 + pow(kappa, 2) ) * ky ) / kappa / pow(1 + ky, 2) : - 2 * CA * ( 1 - kappa ) * ( 1 + kappa - ( 1 - 3 * kappa ) * ky2 ) / kappa / pow(1 - ky2, 2));
  }
  double Cgtmd1gge::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1gge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1gge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1qgo::Cgtmd1qgo(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1qgo::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return 2 * TR / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1ggo::Cgtmd1ggo(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1ggo::LocalPV() const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1)
      return - CA * ( 1 + pow(kappa, 2) ) / pow(kappa, 2) * M_PI;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1polnse::Cgtmd1polnse(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polnse::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) / kappa / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1polnse::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1polqqe::Cgtmd1polqqe(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polqqe::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }
  double Cgtmd1polqqe::Local(double const&) const
  {
    return - CF * zeta2;
  }

  //_________________________________________________________________________________
  Cgtmd1polqge::Cgtmd1polqge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polqge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? 8 * TR * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * TR * ( 3 + ky ) / pow(1 + ky, 2) : 8 * TR * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1polqge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return 2 * TR / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1polqge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return 2 * TR * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1polgqe::Cgtmd1polgqe(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polgqe::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF * ( 1 - y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 4 * CF * ( 1 - kappa ) * y / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1polgge::Cgtmd1polgge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1polgge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * CA * ( 3 + ky ) / pow(1 + ky, 2) : - 8 * CA * ( 1 - kappa ) * y / pow(1 - ky2, 2));
  }
  double Cgtmd1polgge::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1polgge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1polgge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y) / kappa;
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1lingge::Cgtmd1lingge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1lingge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 8 * CA * pow(kappa, 2) * y * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? 2 * CA * kappa * ( 1 - ky ) / pow(1 + ky, 2) : - 8 * CA * ky2 * ( 1 - kappa ) / pow(1 - ky2, 2));
  }
  double Cgtmd1lingge::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1lingge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA * kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1lingge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1linunpgqe::Cgtmd1linunpgqe(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linunpgqe::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF * ( 1 - y ) / y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 4 * CF * ( 1 - kappa ) / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1linunpgge::Cgtmd1linunpgge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linunpgge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 4 * CA * ( 1 - y ) * ( 1 + ky2 ) / y / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - 2 * CA * ( 2 + ky * ( 1 + ky ) ) / y / pow(1 + ky, 2) : - 4 * CA * ( 1 - kappa ) * ( 1 + ky2 ) / pow(1 - ky2, 2));
  }
  double Cgtmd1linunpgge::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1linunpgge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - 2 * CA * kappa / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1linunpgge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - 2 * CA * log(1 - kappa * y);
    else
      return 0;
  }

  //_________________________________________________________________________________
  Cgtmd1linpolgqe::Cgtmd1linpolgqe(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linpolgqe::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 - y ) * kappa / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * ( 1 - kappa ) / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Cgtmd1linpolgge::Cgtmd1linpolgge(double const& xi):
    Expression(1/xi)
  {
  }
  double Cgtmd1linpolgge::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky    = kappa * y;
    const double ky2   = pow(ky, 2);
    if (kappa < 1)
      return (y <= 1 ? - 4 * CA * kappa * ( 1 - y ) / pow(1 - ky2, 2) : 0);
    else
      return (y <= 1 ? - CA * ( 3 + ky ) / pow(1 + ky, 2) : - 4 * CA * ( 1 - kappa ) / pow(1 - ky2, 2));
  }
  double Cgtmd1linpolgge::Local(double const&) const
  {
    return - CA * zeta2;
  }
  double Cgtmd1linpolgge::SingularPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1)
      return - CA / ( 1 - kappa * y );
    else
      return 0;
  }
  double Cgtmd1linpolgge::LocalLogPV(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    if (kappa > 1 && y < 1 && kappa * y < 1)
      return - CA * log(1 - kappa * y) / kappa;
    else
      return 0;
  }
}
