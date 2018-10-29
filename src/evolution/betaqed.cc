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
    // Number of colours
    const int Nc = 3;

    // Sum of the squared electric charges
    const std::vector<double> SumCh2{0., 1./9., 5./9., 2./3., 10./9., 11./9., 5./3.};
    const double coeff = - 4. / 3. * ( Nc * SumCh2[nf] + nl );
    return coeff;
  }

  //_________________________________________________________________________
  double beta1qed(int const& nf, int const& nl)
  {
    // Number of colours
    const int Nc = 3;

    // Sum of the electric charges to the fourth
    const std::vector<double> SumCh4{0., 1./81., 17./81., 18./81., 34./81., 35./81., 51./81.};

    const double coeff = - 16. / 4. * ( Nc * SumCh4[nf] + nl );
    return coeff;
  }
}
