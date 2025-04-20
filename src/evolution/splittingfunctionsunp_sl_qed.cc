//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionsunp_sl_qed.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________________
  P01qedns::P01qedns():
    Expression()
  {
  }
  double P01qedns::Regular(double const& x) const
  {
    return - 2 * ( 1 + x );
  }
  double P01qedns::Singular(double const& x) const
  {
    return 4 / ( 1 - x );
  }
  double P01qedns::Local(double const& x) const
  {
    return 4 * log(1 - x) + 3;
  }

  //_________________________________________________________________________________
  P01qedqgm::P01qedqgm():
    Expression()
  {
  }
  double P01qedqgm::Regular(double const& x) const
  {
    return 4 * ( 1 - 2 * x + 2 * x * x );
  }

  //_________________________________________________________________________________
  P01qedgmq::P01qedgmq():
    Expression()
  {
  }
  double P01qedgmq::Regular(double const& x) const
  {
    return 4 * ( - 1 + 0.5 * x + 1 / x );
  }

  //_________________________________________________________________________________
  P01qedgmgm::P01qedgmgm():
    Expression()
  {
  }
  double P01qedgmgm::Local(double const&) const
  {
    return - 4. / 3;
  }

  //_________________________________________________________________________________
  P02qednsp::P02qednsp(double const& crat2):
    Expression(),
    _crat2(crat2)
  {
    _a2 = - 80 * _crat2 / 9.;
  }
  double P02qednsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1 + x) - Pi2 / 6;
    const double gqq1  =
      + 4 * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x ) - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
      + 4 * _crat2 * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 )
      + 4 * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    return gqq1 - gqq1l;
  }
  double P02qednsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P02qednsp::Local(double const& x) const
  {
    const double p1delta = 3 / 2. + 24 * zeta3 - 2 * Pi2 + _crat2 * ( - 2 / 3. - 8 * Pi2 / 9. );
    return log(1 - x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P02qednsm::P02qednsm(double const& crat2):
    Expression(),
    _crat2(crat2)
  {
    _a2 = - 80 * _crat2 / 9.;
  }
  double P02qednsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqq   = 2 / ( 1 - x ) - 1 - x;
    const double pqqmx = 2 / ( 1 + x ) - 1 + x;
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1 + x) - Pi2 / 6;
    const double gqq1  =
      + 4 * ( ( - 3 * lnx / 2 - 2 * ln1mx * lnx ) * pqq - 5 * ( 1 - x ) - lnx2 * ( 1 + x ) / 2 - lnx * ( 1.5 + 7 * x / 2 ) )
      + 4 * _crat2 * ( ( - 10 / 9. - 2 * lnx / 3 ) * pqq - 4 * ( 1 - x ) / 3 )
      - 4 * ( 2 * pqqmx * S2x + 4 * ( 1 - x ) + 2 * lnx * ( 1 + x ) );
    const double gqq1l = _a2 / ( 1 - x );
    return gqq1 - gqq1l;
  }
  double P02qednsm::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P02qednsm::Local(double const& x) const
  {
    const double p1delta = 3 / 2. + 24 * zeta3 - 2 * Pi2 + _crat2 * ( - 2 / 3. - 8 * Pi2 / 9. );
    return log(1 - x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P02qedps::P02qedps():
    Expression()
  {
  }
  double P02qedps::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return - 8 + 24 * x - 224 / 9. * x * x + 80 / 9. / x
           + ( 4 + 20 * x ) * lnx + 32 / 3. * x * x * lnx
           - ( 4 + 4 * x ) * lnx2;
  }

  //_________________________________________________________________________________
  P02qedqgm::P02qedqgm():
    Expression()
  {
  }
  double P02qedqgm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pqg   = x * x + ( 1 - x ) * ( 1 - x );
    return 4 * ( 4  + 4 * ln1mx + ( 10 - 4 * ( ln1mx - lnx ) + 2 * ( - ln1mx + lnx ) * ( - ln1mx + lnx ) - 2 * Pi2 / 3 ) * pqg
                 - lnx * ( 1 - 4 * x ) - lnx2 * ( 1  - 2 * x ) - 9 * x );
  }

  //_________________________________________________________________________________
  P02qedgmq::P02qedgmq(double const& crat2):
    Expression(),
    _crat2(crat2)
  {
  }
  double P02qedgmq::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1 - x);
    const double pgq   = ( 1 + ( 1 - x ) * ( 1 - x ) ) / x;
    return 4 * ( - 2.5 - ( 3 * ln1mx + ln1mx * ln1mx ) * pgq - lnx2 * ( 1 - x / 2 ) - 7 * x / 2
                 - 2 * ln1mx * x + lnx * ( 2 + 7 * x / 2 ) )
           + 4 * _crat2 * ( - ( 20 / 9. + 4 * ln1mx / 3 ) * pgq - 4 * x / 3 );
  }

  //_________________________________________________________________________________
  P02qedgmgm::P02qedgmgm():
    Expression()
  {
  }
  double P02qedgmgm::Regular(double const& x) const
  {
    const double lnx  = log(x);
    const double lnx2 = lnx * lnx;
    return 4 * ( - 16 + 4 / ( 3 * x ) + 8 * x + ( 20 * x * x ) / 3 - lnx2 * ( 2 + 2 * x ) - lnx * ( 6 + 10 * x ) );
  }
  double P02qedgmgm::Local(double const&) const
  {
    return - 4.;
  }
}
