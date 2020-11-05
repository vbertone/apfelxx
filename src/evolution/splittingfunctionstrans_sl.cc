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
  /**
   * @brief The LO space-like splitting function for tranversely
   * polarised PDFs. Reference
   * https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
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

  /**
   * @brief The NLO space-like splitting function for tranversely
   * polarised PDFs. Reference
   * https://arxiv.org/pdf/hep-ph/9706511v2.pdf.
   */
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
}
