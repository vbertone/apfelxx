//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdsplittingfunctionstrans_sl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  Pgpd0transns::Pgpd0transns(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0transns::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + 2 * kappa * y + kappa * ( 2 - kappa ) * pow(y, 2) ) / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0transns::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0transns::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log(1 - y) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))));
  }
  double Pgpd0transns::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log(1 - y) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0transqq::Pgpd0transqq(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0transqq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 4 * CF / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + 2 * y + pow(kappa, 2) * pow(y, 2) ) / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0transqq::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0transqq::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log(1 - y) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))) );
  }
  double Pgpd0transqq::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log(1 - y) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0transgg::Pgpd0transgg(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd0transgg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    const double k2    = kappa * kappa;
    return 4 * CA * (y <= 1 ? ( k2 * y - 1 ) * ( ky2 + 1 ) / pow(1 - ky2, 2) : 0)
           + 2 * CA * (kappa > 1 ? ( 1 + 2 * y - 4 * ky2 * ( kappa - 1 ) - pow(kappa * y, 3) * ( kappa * y + 4 ) + 2 * ky2 * y ) / pow(1 - ky2, 2) / ( 1 + y ) : 0);
  }
  double Pgpd0transgg::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0transgg::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 4 * CA * log(1 - y) + ( 11 * CA - 4 * _nf * TR ) / 3 - 2 * CA * log(std::abs(1 - pow(kappa, 2)));
  }
  double Pgpd0transgg::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( - (kappa > 1 ? log(1 - y) : 0) );
  }
}
