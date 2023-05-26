//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionstrans_tl.h"
#include "apfel/constants.h"
#include "apfel/specialfunctions.h"

namespace apfel
{
  //_________________________________________________________________________________
  P0Ttransns::P0Ttransns():
    Expression()
  {
  }
  double P0Ttransns::Regular(double const&) const
  {
    return - 4 * CF;
  }
  double P0Ttransns::Singular(double const& x) const
  {
    return 4 * CF / ( 1 - x );
  }
  double P0Ttransns::Local(double const& x) const
  {
    return 4 * CF * log( 1 - x ) + 3 * CF;
  }

  //_________________________________________________________________________________
  P1Ttransnsp::P1Ttransnsp(int const& nf):
    Expression(),
    _nf(nf)
  {
    _a2 = 2 * CA * CF * ( 134. / 9. - 2 * Pi2 / 3 ) - 80 * _nf * CF * TR / 9;
  }
  double P1Ttransnsp::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double func  = 2 * x * lnx / ( 1 - x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * CF * CF * ( 1 - x + ( 3. / 2. + 2 * ln1mx - 2 * lnx ) * func )
      + 2 * CF * CA * ( - 143. / 9. + 2 * Pi2 / 3 + x + ( 11. / 3. + lnx ) * func )
      + 8 * _nf * CF * TR * ( - func + 10. / 3. ) / 3
      + 4 * CF * ( CF - CA / 2 ) * ( - 1 + x - 4 * S2x / ( 1 + x ) );
  }
  double P1Ttransnsp::Singular(double const& x) const
  {
    return _a2 / ( 1 - x );
  }
  double P1Ttransnsp::Local(double const& x) const
  {
    const double p1delta =
      + 4 * CF * CF * ( 3. / 8. - Pi2 / 2 + 6 * zeta3 )
      + 2 * CF * CA * ( 17. / 12. + 11 * Pi2 / 9 - 6 * zeta3 )
      - 8 * _nf * CF * TR * ( 1. / 4. + Pi2 / 3 );
    return log(1-x) * _a2 + p1delta;
  }

  //_________________________________________________________________________________
  P1Ttransnsm::P1Ttransnsm(int const& nf):
    P1Ttransnsp(nf)
  {
  }
  double P1Ttransnsm::Regular(double const& x) const
  {
    const double lnx   = log(x);
    const double lnx2  = lnx * lnx;
    const double ln1mx = log(1-x);
    const double func  = 2 * x * lnx / ( 1 - x );
    const double S2x   = - 2 * dilog(-x) + lnx2 / 2 - 2 * lnx * log(1+x) - Pi2 / 6;
    return
      + 4 * CF * CF * ( 1 - x + ( 3. / 2. + 2 * ln1mx - 2 * lnx ) * func )
      + 2 * CF * CA * ( - 143. / 9. + 2 * Pi2 / 3 + x + ( 11. / 3. + lnx ) * func )
      + 8 * _nf * CF * TR * ( - func + 10. / 3. ) / 3
      - 4 * CF * ( CF - CA / 2 ) * ( - 1 + x - 4 * S2x / ( 1 + x ) );
  }
}
