//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/djet.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double dJetqCone1()
  {
    return CF * ( 7 + 6 * log(2) - 5 * Pi2 / 6 );
  }

  //_________________________________________________________________________
  double dJetqkT1()
  {
    return CF * ( 13 - 3 * Pi2 / 2 );
  }
}
