//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gpdsplittingfunctionspol_sl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  Pgpd0polns::Pgpd0polns(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0polns::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 + y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + ( 1 + kappa ) * y + ( 1 + kappa - pow(kappa, 2) ) * pow(y, 2) ) / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0polns::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0polns::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log(1 - y) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))));
  }
  double Pgpd0polns::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log(1 - y) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0polqq::Pgpd0polqq(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0polqq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( 1 + y ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * ( 1 + y + kappa * y + pow(kappa, 3) * pow(y, 2) ) / kappa / ( 1 + y ) / ( 1 - ky2 ) : 0);
  }
  double Pgpd0polqq::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0polqq::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( 2 * log(1 - y) + 3. / 2. - log(std::abs(1 - pow(kappa, 2))) );
  }
  double Pgpd0polqq::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CF * ( - (kappa > 1 ? log(1 - y) : 0) );
  }

  //_________________________________________________________________________________
  Pgpd0polqg::Pgpd0polqg(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd0polqg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? 4 * _nf * TR * ( 2 * y - 1 - ky2 ) / pow(1 - ky2, 2) : 0)
           + (kappa > 1 ? 4 * _nf * TR * ( - 1 + kappa ) * ( ky2 + 1 ) / kappa / pow(1 - ky2, 2) : 0);
  }

  //_________________________________________________________________________________
  Pgpd0polgq::Pgpd0polgq(double const& xi):
    Expression(1 / xi)
  {
  }
  double Pgpd0polgq::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    return (y <= 1 ? - 2 * CF * ( pow(kappa, 2) * y + y - 2 ) / ( 1 - ky2 ) : 0)
           + (kappa > 1 ? 2 * CF * pow(1 - kappa, 2) / kappa / ( 1 - ky2 ) : 0);
  }

  //_________________________________________________________________________________
  Pgpd0polgg::Pgpd0polgg(int const& nf, double const& xi):
    Expression(1 / xi),
    _nf(nf)
  {
  }
  double Pgpd0polgg::Regular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    const double ky2   = pow(kappa * y, 2);
    const double k2    = kappa * kappa;
    return 4 * CA * (y <= 1 ? ( k2 * y - 1 ) * ( ky2 + 2 * y - 1 ) / pow(1 - ky2, 2) : 0)
           + 2 * CA * (kappa > 1 ? ( - 3 * kappa - kappa * pow(ky2, 2) - 2 * kappa * y * ky2 - k2 * ky2 * ( 1 + y )
                                     + k2 * ( 1 + y ) * ( 1 + y * y ) - 2 * kappa * y + 3 * y + 3 ) / ( kappa * ( 1 + y ) * pow(1 - ky2, 2) ) : 0);
  }
  double Pgpd0polgg::Singular(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( (y <= 1 ? 2 : 0) - (kappa > 1 ? 1 : 0) ) / ( 1 - y );
  }
  double Pgpd0polgg::Local(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 4 * CA * log(1 - y) + ( 11 * CA - 4 * _nf * TR ) / 3 - 2 * CA * log(std::abs(1 - pow(kappa, 2)));
  }
  double Pgpd0polgg::LocalPP(double const& y) const
  {
    const double kappa = 1 / _eta / _extvar;
    return 2 * CA * ( - (kappa > 1 ? log(1 - y) : 0) );
  }
}
