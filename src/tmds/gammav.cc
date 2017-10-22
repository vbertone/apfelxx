//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/gammav.h"
#include "apfel/constants.h"

namespace apfel
{
  //_________________________________________________________________________
  double gammaVq0()
  {
    return - 6 * CF;
  }

  //_________________________________________________________________________
  double gammaVq1(int const& nf)
  {
    const double coeff = CF * CF * ( - 3 + 4 * Pi2 - 48 * zeta3 )
      + CF * CA * ( - 961. / 27. - 11 * Pi2 / 3 + 52 * zeta3 )
      + CF * TR * nf * (  260. / 27. + 4 * Pi2 / 3 );
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVq2(int const& nf)
  {
    const double coeff = pow(CF,3) * (  - 29 - 6 * Pi2 -  16 * Pi2 * Pi2 / 5 - 136 * zeta3 +  32 * Pi2 * zeta3 / 3  + 480 * zeta5 ) 
      + CF * CF * CA * ( - 151. / 2. +  410 * Pi2 / 9 + 494 * Pi2 * Pi2 / 135 - 1688 * zeta3 / 3 - 16 * Pi2 * zeta3 / 3 - 240 * zeta5 ) 
      + CF * CA * CA * ( - 139345. / 1458. - 7163 * Pi2 / 243 - 83 * Pi2 * Pi2 / 45 + 7052 * zeta3 / 9 - 88 * Pi2 * zeta3 / 9 - 272 * zeta5 ) 
      + CF * CF * TR * nf * ( 5906. / 27. - 52 * Pi2 / 9 - 56 * Pi2 * Pi2 / 27 + 1024 * zeta3 / 9 ) 
      + CF * CA * TR * nf * ( - 34636. / 729. + 5188 * Pi2 / 243 + 44 * Pi2 * Pi2 / 45 - 3856 * zeta3 / 27 ) 
      + CF * TR * TR * nf * nf * ( 19336. / 729. - 80 * Pi2 / 27 - 64 * zeta3 / 27 );
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVg0(int const& nf)
  {
    const double coeff = - 22 * CA / 3 + 8 * TR * nf / 3;
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVg1(int const& nf)
  {
    const double coeff = CA * CA * ( - 1384. / 27. + 11 * Pi2 / 9 + 4 * zeta3 )
      + CA * TR * nf * ( 512. / 27. - 4 * Pi2 / 9 ) + 8 * CF * TR * nf;
    return coeff;
  }

  //_________________________________________________________________________
  double gammaVg2(int const& nf)
  {
    const double coeff = 2 * pow(CA,3) * (  - 97186. / 729. + 6109 * Pi2 / 486 - 319  * Pi2 * Pi2 / 270  + 122 * zeta3 / 3 - 20 * Pi2 * zeta3 / 9 - 16 * zeta5 )
      + 2 * CA * CA * TR * nf * ( 30715. / 729. - 1198 * Pi2 / 243 + 82 * Pi2 * Pi2 / 135 + 712 * zeta3 / 27 )
      + 2 * CA * CF * TR * nf * ( 2434. / 27. - 2 * Pi2 / 3 - 8 * Pi2 * Pi2 / 45  - 304 * zeta3 / 9 ) 
      - 4 * CF * CF * TR * nf
      + 2 * CA * TR * TR * nf * nf * ( - 538. / 729. + 40 * Pi2 / 81 - 224 * zeta3 / 27 ) 
      - 88 * CF *  TR * TR * nf * nf / 9;
    return coeff;
  };
}
