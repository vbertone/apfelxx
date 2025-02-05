//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gammam.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double gammam0()
  {
    return 4;
  }

  //_________________________________________________________________________
  double gammam1(int const& nf)
  {
    return 202. / 3. - 20. * nf / 3.;
  }

  //_________________________________________________________________________
  double gammam2(int const& nf)
  {
    return 1249. - ( 2216. / 27. + 160. * zeta3 / 3. ) * nf - 140. / 81. * nf * nf;
  }
}
