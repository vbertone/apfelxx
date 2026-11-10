//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/zeromasscff.h"
#include "apfel/constants.h"

#include <cmath>
#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  ReCFF10ns::ReCFF10ns():
    Expression{}
  {
  }
  double ReCFF10ns::Regular(double const& x)  const
  {
    return - 1 / ( 1 + x );
  }
  double ReCFF10ns::Singular(double const& x) const
  {
    return - 1 / ( 1 - x );
  }
  double ReCFF10ns::LocalPP(double const& c, double const& d) const
  {
    return - ( log(1 - c) - (d > 1 ? log(d - 1) : 0) );
  }

  //_________________________________________________________________________________
  ImCFFL1nsSpaceLike::ImCFFL1nsSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFFL1nsSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 8 * M_PI * CF * y / ( 1 - pow(kappa * y, 2) ) * ( (y < 1 ? 1 : 0) - (kappa > 1 ? 1 / kappa / y : 0) );
  }

  //_________________________________________________________________________________
  ReCFFL1nsSpaceLike::ReCFFL1nsSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFFL1nsSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 8 * CF * y / ( 1 - pow(kappa * y, 2) ) * ( log(std::abs(( 1 - kappa ) / ( 1 + kappa ))) / kappa / y - log(std::abs(( 1 - y ) / ( 1 + y ))) );
  }

  //_________________________________________________________________________________
  ImCFFL1nsTimeLike::ImCFFL1nsTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFFL1nsTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return - 8 * M_PI * CF * y / ( 1 - pow(kappa * y, 2) ) * ( (y < 1 ? 1 : 0) - (kappa > 1 ? 1 / kappa / y : 0) );
  }

  //_________________________________________________________________________________
  ReCFFL1nsTimeLike::ReCFFL1nsTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFFL1nsTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 8 * CF * y / ( 1 - pow(kappa * y, 2) ) * ( log(std::abs(( 1 - kappa ) / ( 1 + kappa ))) / kappa / y - log(std::abs(( 1 - y ) / ( 1 + y ))) );
  }

  //_________________________________________________________________________________
  ImCFFL1gSpaceLike::ImCFFL1gSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFFL1gSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double ts    = (y < 1 ? 1 : 0);
    const double tu    = 0;
    const double tq1   = 0;
    const double tq2   = (kappa > 1 ? 1 : 0);
    return - 16 * M_PI * TR / pow(1 - ky2, 2) * ( 2 * y * ( ( 1 + y ) * tu - ( 1 - y ) * ts )
                                                  - ( ( 1 + ( 2 * kappa + k2 ) * y * y ) * tq1 - ( 1 - ( 2 * kappa - k2 ) * y * y ) * tq2 ) / kappa );
  }

  //_________________________________________________________________________________
  ReCFFL1gSpaceLike::ReCFFL1gSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFFL1gSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    return 16 * TR / pow(1 - ky2, 2) * ( 2 * y * ( ( 1 + y ) * log(( 1 + y ) / y) - ( 1 - y ) * log(std::abs(1 - y) / y) )
                                         - ( ( 1 + ( 2 * kappa + k2 ) * y * y ) * log(1 + kappa) - ( 1 - ( 2 * kappa - k2 ) * y * y ) * log(std::abs(1 - kappa)) ) / kappa );
  }

  //_________________________________________________________________________________
  ImCFFL1gTimeLike::ImCFFL1gTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFFL1gTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    // Interpolate linearly across the double pole at y = 1 / kappa,
    // where the cancellation is numerically unstable
    const double eta = std::min(0.1, 2e-6 / kappa / kappa);
    if (std::abs(1 - kappa * y) < eta && std::abs(1 - kappa) > 4 * eta)
      {
        const double ym = ( 1 - 2 * eta ) / kappa;
        const double yp = ( 1 + 2 * eta ) / kappa;
        const double fm = Regular(ym);
        return fm + ( Regular(yp) - fm ) * ( y - ym ) / ( yp - ym );
      }
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double ts    = (y > 1 ? 1 : 0);
    const double tu    = 1;
    const double tq1   = 1;
    const double tq2   = (kappa < 1 ? 1 : 0);
    return - 16 * M_PI * TR / pow(1 - ky2, 2) * ( 2 * y * ( ( 1 + y ) * tu - ( 1 - y ) * ts )
                                                  - ( ( 1 + ( 2 * kappa + k2 ) * y * y ) * tq1 - ( 1 - ( 2 * kappa - k2 ) * y * y ) * tq2 ) / kappa );
  }

  //_________________________________________________________________________________
  ReCFFL1gTimeLike::ReCFFL1gTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFFL1gTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    // Interpolate linearly across the double pole at y = 1 / kappa,
    // where the cancellation is numerically unstable
    const double eta = std::min(0.1, 2e-6 / kappa / kappa);
    if (std::abs(1 - kappa * y) < eta && std::abs(1 - kappa) > 4 * eta)
      {
        const double ym = ( 1 - 2 * eta ) / kappa;
        const double yp = ( 1 + 2 * eta ) / kappa;
        const double fm = Regular(ym);
        return fm + ( Regular(yp) - fm ) * ( y - ym ) / ( yp - ym );
      }
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    return 16 * TR / pow(1 - ky2, 2) * ( 2 * y * ( ( 1 + y ) * log(( 1 + y ) / y) - ( 1 - y ) * log(std::abs(1 - y) / y) )
                                         - ( ( 1 + ( 2 * kappa + k2 ) * y * y ) * log(1 + kappa) - ( 1 - ( 2 * kappa - k2 ) * y * y ) * log(std::abs(1 - kappa)) ) / kappa );
  }

  //_________________________________________________________________________________
  ImCFF11nsSpaceLike::ImCFF11nsSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFF11nsSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    return 2 * M_PI * CF * ( (y < 1 ? ( 3 - ( 1 + y ) * log(1 - y) ) / ( 1 - ky2 ) - ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 - ky2 ) * log(y) / ( 1 - y ) : 0)
                             + (kappa > 1 ? ( 2 * ( 1 + y + kappa * y + kappa * k2 * y * y ) * log(kappa - 1) - 3 * kappa * ( 1 + 2 * y + ky2 ) ) / ( 2 * kappa * ( 1 + y ) * ( 1 - ky2 ) ) : 0) );
  }
  double ImCFF11nsSpaceLike::Singular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 2 * M_PI * CF * ( (y < 1 ? 2 * log(1 - y) - 3. / 2. : 0) + (kappa > 1 ? 3. / 2. - log(kappa - 1) : 0) ) / ( 1 - y );
  }
  double ImCFF11nsSpaceLike::Local(double const& c) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 2 * M_PI * CF * ( pow(log(1 - c), 2) - 3 * log(1 - c) / 2
                             - ( 9. / 2. + M_PI * M_PI / 3 + pow(log(1 + kappa), 2) / 2 + pow(log(std::abs(1 - kappa)), 2) / 2
                                 - 3 * log(std::abs(1 - kappa * kappa)) / 2 - (kappa > 1 ? M_PI * M_PI / 2 : 0) ) );
  }
  double ImCFF11nsSpaceLike::LocalPP(double const& c, double const& d) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 2 * M_PI * CF * (kappa > 1 ? ( 3. / 2. - log(kappa - 1) ) * ( log(1 - c) - (d > 1 ? log(d - 1) : 0) ) : 0);
  }

  //_________________________________________________________________________________
  ReCFF11nsSpaceLike::ReCFF11nsSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFF11nsSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double pi2   = M_PI * M_PI;
    const double l1y   = log(std::abs(1 - y));
    const double ly    = log(y);
    const double lpy   = log(( 1 + y ) / y);
    const double lpk   = log(1 + kappa);
    const double lmk   = log(std::abs(1 - kappa));
    return CF / 2 * ( ( 2 * ( 1 + y ) * l1y * l1y - 12 * l1y - 8 * y * log(1 + y) - (y < 1 ? 2 * pi2 * ( 1 + y ) : 0) ) / ( 1 - ky2 )
                      + 2 * ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 - ky2 ) * ( 2 * ly * l1y - ly * ly ) / ( 1 - y )
                      + 2 * ( 3 - 2 * y - ( 4 - 3 * k2 ) * y * y ) / ( 1 - ky2 ) * ly / ( 1 - y )
                      - 2 * ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 + y ) / ( 1 - ky2 ) * lpy * lpy
                      - 2 * ( 3 + 2 * y - ( 4 - 3 * k2 ) * y * y ) / ( 1 + y ) / ( 1 - ky2 ) * lpy
                      + 2 / kappa / ( 1 + y ) / ( 1 - ky2 ) * ( ( 1 + y - kappa * y - kappa * k2 * y * y ) * lpk * lpk
                                                                - ( 1 + y + kappa * y + kappa * k2 * y * y ) * ( lmk * lmk - (kappa > 1 ? pi2 : 0) ) )
                      + 6 * ( 1 + 2 * y + ky2 ) / ( 1 + y ) / ( 1 - ky2 ) * log(std::abs(1 - k2))
                      + 18 / ( 1 + y ) );
  }
  double ReCFF11nsSpaceLike::Singular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double pi2   = M_PI * M_PI;
    const double l1y   = log(std::abs(1 - y));
    return CF / 2 * ( - 4 * l1y * l1y + 6 * l1y + (y < 1 ? 4 * pi2 : 0)
                      + 2 * ( 9 - 3 * log(std::abs(1 - kappa * kappa)) + pow(log(1 + kappa), 2) + pow(log(std::abs(1 - kappa)), 2) - (kappa > 1 ? pi2 : 0) ) ) / ( 1 - y );
  }
  double ReCFF11nsSpaceLike::Local(double const& c) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double pi2   = M_PI * M_PI;
    return CF / 2 * ( 4 * pi2 * log(1 - c) - pi2 * ( 3 + (kappa > 1 ? 4 * log(kappa - 1) - 6 : 0) ) );
  }
  double ReCFF11nsSpaceLike::LocalPP(double const& c, double const& d) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double lc    = log(1 - c);
    const double ld    = (d > 1 ? log(d - 1) : 0);
    const double b     = 9 - 3 * log(std::abs(1 - kappa * kappa)) + pow(log(1 + kappa), 2) + pow(log(std::abs(1 - kappa)), 2) - (kappa > 1 ? M_PI * M_PI : 0);
    return CF / 2 * ( - 4 * ( pow(lc, 3) - pow(ld, 3) ) / 3 + 3 * ( pow(lc, 2) - pow(ld, 2) ) + 2 * b * ( lc - ld ) );
  }

  //_________________________________________________________________________________
  ImCFF11nsTimeLike::ImCFF11nsTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFF11nsTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double lpk   = log(1 + kappa);
    return 2 * M_PI * CF * ( (y > 1 ? ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 - ky2 ) * log(y) / ( y - 1 ) - ( 1 + y ) * log(y - 1) / ( 1 - ky2 ) : 0)
                             - (y < 1 ? 3 / ( 1 - ky2 ) : 0)
                             + ( 2 + 3 * kappa - (kappa > 1 ? 2 : 0) ) / kappa / ( 1 - ky2 )
                             + ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 + y ) / ( 1 - ky2 ) * log(( 1 + y ) / y)
                             + (kappa < 1 ? ( 2 * ( 1 + y + kappa * y + kappa * k2 * y * y ) * log(1 - kappa)
                                              - 4 - 4 * y - 3 * kappa - 6 * kappa * y - 3 * kappa * k2 * y * y ) / ( 2 * kappa * ( 1 + y ) * ( 1 - ky2 ) ) : 0)
                             - ( 1 + y - kappa * y - kappa * k2 * y * y ) * lpk / kappa / ( 1 + y ) / ( 1 - ky2 ) );
  }
  double ImCFF11nsTimeLike::Singular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 2 * M_PI * CF * ( 2 * log(std::abs(1 - y)) + (y < 1 ? - 2 * log(1 - y) + 3. / 2. : 0)
                             + (kappa < 1 ? 3. / 2. - log(1 - kappa) : 0) - log(1 + kappa) ) / ( 1 - y );
  }
  double ImCFF11nsTimeLike::Local(double const& c) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    return 2 * M_PI * CF * ( - pow(log(1 - c), 2) + 3 * log(1 - c) / 2
                             + 9. / 2. - M_PI * M_PI / 6 + pow(log(1 + kappa), 2) / 2 + pow(log(std::abs(1 - kappa)), 2) / 2
                             - 3 * log(std::abs(1 - kappa * kappa)) / 2 - (kappa < 1 ? M_PI * M_PI / 2 : 0) );
  }
  double ImCFF11nsTimeLike::LocalPP(double const& c, double const& d) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double lc    = log(1 - c);
    const double ld    = (d > 1 ? log(d - 1) : 0);
    const double a     = (kappa < 1 ? 3. / 2. - log(1 - kappa) : 0) - log(1 + kappa);
    return 2 * M_PI * CF * ( pow(lc, 2) - pow(ld, 2) + a * ( lc - ld ) );
  }

  //_________________________________________________________________________________
  ReCFF11nsTimeLike::ReCFF11nsTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFF11nsTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double pi2   = M_PI * M_PI;
    const double l1y   = log(std::abs(1 - y));
    const double ly    = log(y);
    const double lpy   = log(( 1 + y ) / y);
    const double lpk   = log(1 + kappa);
    const double lmk   = log(std::abs(1 - kappa));
    return CF / 2 * ( ( 2 * ( 1 + y ) * l1y * l1y - 12 * l1y - 8 * y * log(1 + y) + (y < 1 ? 2 * pi2 * ( 1 + y ) : 0) ) / ( 1 - ky2 )
                      + 2 * ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 - ky2 ) * ( 2 * ly * l1y - ly * ly ) / ( 1 - y )
                      + 2 * ( 3 - 2 * y - ( 4 - 3 * k2 ) * y * y ) / ( 1 - ky2 ) * ly / ( 1 - y )
                      - 2 * ( 1 + ( 1 - 2 * k2 ) * y * y ) / ( 1 + y ) / ( 1 - ky2 ) * lpy * lpy
                      - 2 * ( 3 + 2 * y - ( 4 - 3 * k2 ) * y * y ) / ( 1 + y ) / ( 1 - ky2 ) * lpy
                      + 2 / kappa / ( 1 + y ) / ( 1 - ky2 ) * ( ( 1 + y - kappa * y - kappa * k2 * y * y ) * lpk * lpk
                                                                - ( 1 + y + kappa * y + kappa * k2 * y * y ) * ( lmk * lmk + (kappa > 1 ? pi2 : 0) ) )
                      + 6 * ( 1 + 2 * y + ky2 ) / ( 1 + y ) / ( 1 - ky2 ) * log(std::abs(1 - k2))
                      + 18 / ( 1 + y ) );
  }
  double ReCFF11nsTimeLike::Singular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double pi2   = M_PI * M_PI;
    const double l1y   = log(std::abs(1 - y));
    return CF / 2 * ( - 4 * l1y * l1y + 6 * l1y - (y < 1 ? 4 * pi2 : 0)
                      + 2 * ( 9 - 3 * log(std::abs(1 - kappa * kappa)) + pow(log(1 + kappa), 2) + pow(log(std::abs(1 - kappa)), 2) + (kappa > 1 ? pi2 : 0) ) ) / ( 1 - y );
  }
  double ReCFF11nsTimeLike::Local(double const& c) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double pi2   = M_PI * M_PI;
    return CF / 2 * ( - 4 * pi2 * log(1 - c) + pi2 * ( 4 * log(1 + kappa) - 3 + (kappa < 1 ? 4 * log(1 - kappa) - 6 : 0) ) );
  }
  double ReCFF11nsTimeLike::LocalPP(double const& c, double const& d) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double lc    = log(1 - c);
    const double ld    = (d > 1 ? log(d - 1) : 0);
    const double b     = 9 - 3 * log(std::abs(1 - kappa * kappa)) + pow(log(1 + kappa), 2) + pow(log(std::abs(1 - kappa)), 2) + (kappa > 1 ? M_PI * M_PI : 0);
    return CF / 2 * ( - 4 * ( pow(lc, 3) - pow(ld, 3) ) / 3 + 3 * ( pow(lc, 2) - pow(ld, 2) ) + 2 * b * ( lc - ld ) );
  }

  //_________________________________________________________________________________
  ImCFF11gSpaceLike::ImCFF11gSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFF11gSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double ts    = (y < 1 ? 1 : 0);
    const double tu    = 0;
    const double tq1   = 0;
    const double tq2   = (kappa > 1 ? 1 : 0);
    const double lsy   = log(std::abs(1 - y) / y);
    const double lpy   = log(( 1 + y ) / y);
    return 4 * M_PI * TR / pow(1 - ky2, 2) * ( ts * ( ( y * y + pow(1 - y, 2) - ky2 ) * lsy - 1 + 4 * y - ( 4 - k2 ) * y * y )
                                               + tu * ( ( pow(1 + y, 2) + y * y - ky2 ) * lpy - 1 - 4 * y - ( 4 - k2 ) * y * y )
                                               - tq1 / kappa * ( ( 1 + kappa ) * ( 1 + ( 2 * kappa - k2 ) * y * y ) * log(1 + kappa) - 1 - kappa - ( 4 * kappa + 3 * k2 - kappa * k2 ) * y * y )
                                               + tq2 / kappa * ( ( 1 - kappa ) * ( 1 - ( 2 * kappa + k2 ) * y * y ) * log(std::abs(1 - kappa)) - 1 + kappa + ( 4 * kappa - 3 * k2 - kappa * k2 ) * y * y ) );
  }

  //_________________________________________________________________________________
  ReCFF11gSpaceLike::ReCFF11gSpaceLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFF11gSpaceLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double ts    = (y < 1 ? 1 : 0);
    const double tu    = 0;
    const double tq1   = 0;
    const double tq2   = (kappa > 1 ? 1 : 0);
    const double pi2   = M_PI * M_PI;
    const double lsy   = log(std::abs(1 - y) / y);
    const double lpy   = log(( 1 + y ) / y);
    const double lpk   = log(1 + kappa);
    const double lmk   = log(std::abs(1 - kappa));
    return 2 * TR / pow(1 - ky2, 2) * ( - ( y * y + pow(1 - y, 2) - ky2 ) * ( lsy * lsy - pi2 * ts )
                                        + 2 * ( 1 - 4 * y + ( 4 - k2 ) * y * y ) * lsy
                                        - ( pow(1 + y, 2) + y * y - ky2 ) * ( lpy * lpy - pi2 * tu )
                                        + 2 * ( 1 + 4 * y + ( 4 - k2 ) * y * y ) * lpy
                                        + ( ( 1 + kappa ) * ( 1 + ( 2 * kappa - k2 ) * y * y ) * ( lpk * lpk - pi2 * tq1 )
                                            - 2 * ( 1 + kappa + ( 4 * kappa + 3 * k2 - kappa * k2 ) * y * y ) * lpk
                                            - ( 1 - kappa ) * ( 1 - ( 2 * kappa + k2 ) * y * y ) * ( lmk * lmk - pi2 * tq2 )
                                            + 2 * ( 1 - kappa - ( 4 * kappa - 3 * k2 - kappa * k2 ) * y * y ) * lmk ) / kappa );
  }

  //_________________________________________________________________________________
  ImCFF11gTimeLike::ImCFF11gTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ImCFF11gTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    // Interpolate linearly across the double pole at y = 1 / kappa,
    // where the cancellation is numerically unstable
    const double eta = std::min(0.1, 2e-6 / kappa / kappa);
    if (std::abs(1 - kappa * y) < eta && std::abs(1 - kappa) > 4 * eta)
      {
        const double ym = ( 1 - 2 * eta ) / kappa;
        const double yp = ( 1 + 2 * eta ) / kappa;
        const double fm = Regular(ym);
        return fm + ( Regular(yp) - fm ) * ( y - ym ) / ( yp - ym );
      }
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double ts    = (y > 1 ? 1 : 0);
    const double tu    = 1;
    const double tq1   = 1;
    const double tq2   = (kappa < 1 ? 1 : 0);
    const double lsy   = log(std::abs(1 - y) / y);
    const double lpy   = log(( 1 + y ) / y);
    return 4 * M_PI * TR / pow(1 - ky2, 2) * ( ts * ( ( y * y + pow(1 - y, 2) - ky2 ) * lsy - 1 + 4 * y - ( 4 - k2 ) * y * y )
                                               + tu * ( ( pow(1 + y, 2) + y * y - ky2 ) * lpy - 1 - 4 * y - ( 4 - k2 ) * y * y )
                                               - tq1 / kappa * ( ( 1 + kappa ) * ( 1 + ( 2 * kappa - k2 ) * y * y ) * log(1 + kappa) - 1 - kappa - ( 4 * kappa + 3 * k2 - kappa * k2 ) * y * y )
                                               + tq2 / kappa * ( ( 1 - kappa ) * ( 1 - ( 2 * kappa + k2 ) * y * y ) * log(std::abs(1 - kappa)) - 1 + kappa + ( 4 * kappa - 3 * k2 - kappa * k2 ) * y * y ) );
  }

  //_________________________________________________________________________________
  ReCFF11gTimeLike::ReCFF11gTimeLike(double const& xi):
    Expression(1 / xi)
  {
  }
  double ReCFF11gTimeLike::Regular(double const& y) const
  {
    const double kappa = std::max(1 / _eta / _extvar, eps8);
    // Interpolate linearly across the double pole at y = 1 / kappa,
    // where the cancellation is numerically unstable
    const double eta = std::min(0.1, 2e-6 / kappa / kappa);
    if (std::abs(1 - kappa * y) < eta && std::abs(1 - kappa) > 4 * eta)
      {
        const double ym = ( 1 - 2 * eta ) / kappa;
        const double yp = ( 1 + 2 * eta ) / kappa;
        const double fm = Regular(ym);
        return fm + ( Regular(yp) - fm ) * ( y - ym ) / ( yp - ym );
      }
    const double k2    = kappa * kappa;
    const double ky2   = k2 * y * y;
    const double ts    = (y > 1 ? 1 : 0);
    const double tu    = 1;
    const double tq1   = 1;
    const double tq2   = (kappa < 1 ? 1 : 0);
    const double pi2   = M_PI * M_PI;
    const double lsy   = log(std::abs(1 - y) / y);
    const double lpy   = log(( 1 + y ) / y);
    const double lpk   = log(1 + kappa);
    const double lmk   = log(std::abs(1 - kappa));
    return 2 * TR / pow(1 - ky2, 2) * ( - ( y * y + pow(1 - y, 2) - ky2 ) * ( lsy * lsy - pi2 * ts )
                                        + 2 * ( 1 - 4 * y + ( 4 - k2 ) * y * y ) * lsy
                                        - ( pow(1 + y, 2) + y * y - ky2 ) * ( lpy * lpy - pi2 * tu )
                                        + 2 * ( 1 + 4 * y + ( 4 - k2 ) * y * y ) * lpy
                                        + ( ( 1 + kappa ) * ( 1 + ( 2 * kappa - k2 ) * y * y ) * ( lpk * lpk - pi2 * tq1 )
                                            - 2 * ( 1 + kappa + ( 4 * kappa + 3 * k2 - kappa * k2 ) * y * y ) * lpk
                                            - ( 1 - kappa ) * ( 1 - ( 2 * kappa + k2 ) * y * y ) * ( lmk * lmk - pi2 * tq2 )
                                            + 2 * ( 1 - kappa - ( 4 * kappa - 3 * k2 - kappa * k2 ) * y * y ) * lmk ) / kappa );
  }
}
