//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/betaqcdqed.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double beta1qcdqed(int const& nf)
  {
    return - 4 * TR * SumCh2[nf];
  }

  //_________________________________________________________________________
  double beta1qedqcd(int const& nf)
  {
    return - 4 * CF * NC * SumCh2[nf];
  }
}
