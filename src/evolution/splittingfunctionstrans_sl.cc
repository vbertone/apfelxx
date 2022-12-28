//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionstrans_sl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

namespace apfel
{
  //_________________________________________________________________________________
  P0transns::P0transns():
    Expression()
  {
  }
  double P0transns::Regular(double const&) const
  {
    return - 4 * CF;
  }
  double P0transns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double P0transns::Local(double const& x) const
  {
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }

  //_________________________________________________________________________________
  P0transgg::P0transgg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0transgg::Regular(double const&) const
  {
    return - 4 * CA;
  }
  double P0transgg::Singular(double const& x) const
  {
    return 4 * CA / ( 1 - x );
  }
  double P0transgg::Local(double const& x) const
  {
    return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
  }

  //_________________________________________________________________________________
  P1transnsp::P1transnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = 2 * CA * CF * ( 134. / 9. - 2 * Pi2 / 3 ) - 80 * _nf * CF * TR / 9;
  }
  double P1transnsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double func  = 2 * x * lnx / ( 1 - x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * CF * CF * ( 1 - x - ( 3. / 2. + 2 * ln1mx ) * func )
      + 2 * CF * CA * ( - 143. / 9. + 2 * Pi2 / 3 + x + ( 11. / 3. + lnx ) * func )
      + 8 * _nf * CF * TR * ( - func + 10. / 3. ) / 3
      + 4 * CF * ( CF - CA / 2 ) * ( - 1 + x - 4 * S2x / ( 1 + x ) );
  }
  double P1transnsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1transnsp::Local(double const& x) const
  {
    const double p1delta =
      + 4 * CF * CF * ( 3. / 8. - Pi2 / 2 + 6 * zeta3 )
      + 2 * CF * CA * ( 17. / 12. + 11 * Pi2 / 9 - 6 * zeta3 )
      - 8 * _nf * CF * TR * ( 1. / 4. + Pi2 / 3 );
    return log(1-x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P1transnsm::P1transnsm(int const& nf):
    P1transnsp(nf)
  {
  }
  double P1transnsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double func  = 2 * x * lnx / ( 1 - x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * CF * CF * ( 1 - x - ( 3. / 2. + 2 * ln1mx ) * func )
      + 2 * CF * CA * ( - 143. / 9. + 2 * Pi2 / 3 + x + ( 11. / 3. + lnx ) * func )
      + 8 * _nf * CF * TR * ( - func + 10. / 3. ) / 3
      - 4 * CF * ( CF - CA / 2 ) * ( - 1 + x - 4 * S2x / ( 1 + x ) );
  }

  //_________________________________________________________________________________
  P1transgg::P1transgg(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = 2 * CA * CA * ( 134. / 9. - 2 * Pi2 / 3 ) - 80 * _nf * CA * TR / 9;
  }
  double P1transgg::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double ln1mx = log(1-x);
    const double x3    = pow(x, 3);
    const double S2x   = - 2 * dilog(-x) + lnx * lnx / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * CA * CA * ( x * lnx * ( lnx - 4 * ln1mx ) / ( 1 - x ) - 67. / 9. + Pi2 / 3 + ( 1 - x3 ) / 6 / x + S2x * ( 2 / ( 1 + x ) - 2 ) )
      + 4 * _nf * TR * CA * ( 20. / 9. + ( 1 - x3 ) / 3 / x )
      - 8 * _nf * TR * CF * ( 1 - x3 ) / 3 / x;
  }
  double P1transgg::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1transgg::Local(double const& x) const
  {
    const double p1delta =
      + 4 * CA * CA * ( 8. / 3. + 3 * zeta3 )
      - 16 * _nf * TR * CA / 3
      - 4 * _nf * TR * CF;
    return log(1-x) * _a2 + p1delta;
  }
}
