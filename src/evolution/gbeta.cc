//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/gbeta.h"
#include "apfel/betaqcd.h"

#include <math.h>

namespace apfel
{
  //_________________________________________________________________________________
  double g1beta(double const& lambda)
  {
    return 1 / ( 1 - lambda );
  }

  //_________________________________________________________________________________
  double g2beta(int const& nf, double const& kappa, double const& lambda)
  {
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1  = - 2 * beta1qcd(nf) / bt0;
    return - ( b1 * log(1 - lambda) + bt0 * log(kappa) ) / pow(1 - lambda, 2);
  }

  //_________________________________________________________________________________
  double g3beta(int const& nf, double const& kappa, double const& lambda)
  {
    const double lnk = log(kappa);
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1  = - 2 * beta1qcd(nf) / bt0;
    const double b2  = - 2 * beta2qcd(nf) / bt0;
    const double ln1ml = log(1 - lambda);
    return ( b2 * lambda - pow(b1, 2) * ( lambda + ln1ml - pow(ln1ml, 2) )
             + bt0 * b1 * ( 2 * ln1ml - 1 ) * lnk + pow(bt0 * lnk, 2) ) / pow(1 - lambda, 3);
  }

  //_________________________________________________________________________________
  double g4beta(int const& nf, double const& kappa, double const& lambda)
  {
    const double lnk = log(kappa);
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1  = - 2 * beta1qcd(nf) / bt0;
    const double b2  = - 2 * beta2qcd(nf) / bt0;
    const double b3  = - 2 * beta3qcd(nf) / bt0;
    const double ln1ml = log(1 - lambda);
    return ( ( b3 - b2 * b1 ) * lambda - ( pow(b1, 3) - 2 * b2 * b1 + b3 ) * pow(lambda, 2) / 2
             + ( 2 * b1 * ( pow(b1, 2) - b2 ) * lambda - b2 * b1 ) * ln1ml
             + pow(b1, 3 ) * ( 5. / 2. - ln1ml ) * pow(ln1ml, 2)
             + ( 2 * ( pow(b1, 2) - b2 ) * lambda + pow(b1, 2) * ( 5 - 3 * ln1ml ) * ln1ml - b2 ) * bt0 * lnk
             + b1 * ( - 3 * ln1ml + 5. / 2. ) * pow(bt0 * lnk, 2) - pow(bt0 * lnk, 3)
           ) / pow(1 - lambda, 4);
  }
}
