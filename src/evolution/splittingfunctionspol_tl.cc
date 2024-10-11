//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionspol_tl.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________________
  P0Tpolns::P0Tpolns():
    P0Tns()
  {
  }

  //_________________________________________________________________________________
  P0Tpolqg::P0Tpolqg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0Tpolqg::Regular(double const& x) const
  {
    return 4 * _nf * CF * ( 2 - x );
  }

  //_________________________________________________________________________________
  P0Tpolgq::P0Tpolgq():
    Expression()
  {
  }
  double P0Tpolgq::Regular(double const& x) const
  {
    return 2 * TR * ( 2 * x - 1 );
  }

  //_________________________________________________________________________________
  P0Tpolgg::P0Tpolgg(int const& nf):
    Expression(),
    _nf(nf)
  {
  }
  double P0Tpolgg::Regular(double const& x) const
  {
    return 4 * CA * ( - 2 * x + 1 );
  }
  double P0Tpolgg::Singular(double const& x) const
  {
    return 4 * CA / ( 1 - x );
  }
  double P0Tpolgg::Local(double const& x) const
  {
    return 4 * CA * log( 1 - x ) - 2 / 3. * _nf + 11 / 3. * CA;
  }
}
