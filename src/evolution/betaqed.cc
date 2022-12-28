//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/betaqed.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double beta0qed(int const& nf, int const& nl)
  {
    return - 4. / 3. * ( NC * SumCh2[nf] + nl );
  }

  //_________________________________________________________________________
  double beta1qed(int const& nf, int const& nl)
  {
    return - 16. / 4. * ( NC * SumCh4[nf] + nl );
  }
}
