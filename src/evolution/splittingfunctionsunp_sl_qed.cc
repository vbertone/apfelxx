//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/splittingfunctionsunp_sl_qed.h"

namespace apfel
{
  //_________________________________________________________________________________
  P0qedns::P0qedns():
    Expression()
  {
  }
  double P0qedns::Regular(double const& x) const
  {
    return - 2 * ( 1 + x );
  }
  double P0qedns::Singular(double const& x) const
  {
    return 4 / ( 1 - x );
  }
  double P0qedns::Local(double const& x) const
  {
    return 4 * log(1 - x) + 3;
  }

  //_________________________________________________________________________________
  P0qedqgm::P0qedqgm():
    Expression()
  {
  }
  double P0qedqgm::Regular(double const& x) const
  {
    return 4 * ( 1 - 2 * x + 2 * x * x );
  }

  //_________________________________________________________________________________
  P0qedgmq::P0qedgmq():
    Expression()
  {
  }
  double P0qedgmq::Regular(double const& x) const
  {
    return 4 * ( - 1 + 0.5 * x + 1 / x );
  }

  //_________________________________________________________________________________
  P0qedgmgm::P0qedgmgm():
    Expression()
  {
  }
  double P0qedgmgm::Local(double const&) const
  {
    return - 4. / 3.;
  }
}
