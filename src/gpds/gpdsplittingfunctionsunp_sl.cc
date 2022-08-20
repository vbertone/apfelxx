//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdsplittingfunctionsunp_sl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  Pgpd0ns::Pgpd0ns(double const& xi):
    Expression(1/xi)
  {
  }
  double Pgpd0ns::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 + y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + ( 1 + kappa ) * y + ( 1 + kappa - pow(kappa, 2) ) * pow(y, 2) ) / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0ns::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0ns::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log( 1 - y ) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))) );
  }
  double Pgpd0ns::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log( 1 - y ) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0qq::Pgpd0qq(double const& xi):
    Expression(1/xi)
  {
  }
  double Pgpd0qq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 + y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + y + kappa * y + pow(kappa, 3) * pow(y, 2) ) / kappa / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0qq::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0qq::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log( 1 - y ) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))) );
  }
  double Pgpd0qq::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log( 1 - y ) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0qg::Pgpd0qg(int const& nf, double const& xi):
    Expression(1/xi),
    _nf(nf)
  {
  }
  double Pgpd0qg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 4 * _nf * TR * ( pow(y, 2) + pow(1 - y, 2) - ky2 ) / pow(1 - ky2, 2) : 0)
           + (kappa > 1 ? 4 * _nf * TR * ( 1 - kappa ) * ( 1 - kappa * ( kappa + 2 ) * pow(y, 2) ) / kappa / pow(1 - ky2, 2) : 0);
  }

  //_________________________________________________________________________________
  Pgpd0gq::Pgpd0gq(double const& xi):
    Expression(1/xi)
  {
  }
  double Pgpd0gq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 2 * CF * ( 1 + pow(1 - y, 2) - ky2 ) / y / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? - 2 * CF * pow(1 - kappa, 2) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Pgpd0gg::Pgpd0gg(int const& nf, double const& xi):
    Expression(1/xi),
    _nf(nf)
  {
  }
  double Pgpd0gg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    const double k2    = kappa * kappa;
    return 4 * CA * (y <= 1 ? - ( 1 + k2 * y ) / ( 1 - ky2 ) + ( ( 1 - y ) / y + y * ( 1 - y ) ) / pow(1 - ky2, 2) : 0)
           + 2 * CA * (kappa > 1 ? 2 * ( 1 - kappa ) * ( 1 + pow(y, 2) ) / pow(1 - ky2, 2) + k2 * ( 1 + y ) / ( 1 - ky2 )
                       + ( 1 - k2 ) * ( 2 - 1 / kappa - 1 / ( 1 + y ) ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0gg::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0gg::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 4 * CA * log( 1 - y ) + ( 11 * CA - 4 * _nf * TR ) / 3 - 2 * CA * log(std::abs(1 - pow(kappa, 2)));
  }
  double Pgpd0gg::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( - (kappa > 1 ? log( 1 - y ) : 0) );
  }
}
