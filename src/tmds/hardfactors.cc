//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/hardfactors.h"
#include "apfel/constants.h"

namespace apfel {
  //_____________________________________________________________________________
  double HardFactorDY(int const& pt, double const& Alphas, int const& nf, double const& kappa)
    {
      // Compute log anf its powers.
      const double lQ2  = 2 * log(kappa);
      const double lQ22 = lQ2 * lQ2;
      const double lQ23 = lQ22 * lQ2;
      const double lQ24 = lQ23 * lQ2;

      // Compute coupling and its powers.
      const double as  = Alphas / FourPi;
      const double as2 = as * as;

      // Now compute hard factor according to the perturbative order.
      double hfct = 1;
      if (pt > 0)
	hfct += 2 * as * CF * ( - lQ22 - 3 * lQ2 - 8 + 7 * Pi2 / 6 );
      if (pt > 1)
	hfct += 2 * as2 * CF *
	  ( CF * ( lQ24 + 6 * lQ23 + ( 25. - 7 * Pi2 / 3 ) * lQ22 + ( 93. / 2. - 5 * Pi2 - 24 * zeta3 ) * lQ2
		   + 511. / 8. - 83 * Pi2 / 6 - 30 * zeta3 + 67 * Pi2 * Pi2 / 60 ) +
	    CA * ( - 11 * lQ23 / 9 + ( - 233. / 18. + Pi2 / 3 ) * lQ22 + ( - 2545. / 54. + 22 * Pi2 / 9 + 26 * zeta3 ) * lQ2
		   - 51157. / 648. + 1061 * Pi2 / 108 + 313 * zeta3 / 9 - 4 * Pi2 * Pi2 / 45 ) +
	    TR * nf * ( 4 * lQ23 / 9 + 38 * lQ22 / 9 + ( 418. / 27. - 8 * Pi2 / 9 ) * lQ2
			+ 4085. / 162. - 91 * Pi2 / 27 + 4 * zeta3 / 9 ) );

      return hfct;
    };

  //_____________________________________________________________________________
  double HardFactorSIDIS(int const& pt, double const& Alphas, int const& nf, double const& kappa)
    {
      // Compute log anf its powers.
      const double lQ2  = 2 * log(kappa);
      const double lQ22 = lQ2 * lQ2;
      const double lQ23 = lQ22 * lQ2;
      const double lQ24 = lQ23 * lQ2;

      // Compute coupling and its powers.
      const double as  = Alphas / FourPi;
      const double as2 = as * as;

      // Now compute hard factor according to the perturbative order.
      double hfct = 1;
      if (pt > 0)
	hfct += 2 * as * CF * ( - lQ22 - 3 * lQ2 - 8 + zeta2 );
      if (pt > 1)
	hfct += 2 * as2 * CF *
	  ( CF * ( lQ24 + 6 * lQ23 + ( 25. - 2 * zeta2 ) * lQ22 + ( 93. / 2. + 6 * zeta2 - 24 * zeta3 ) * lQ2
		   + 511. / 8. + 13 * zeta2 - 30 * zeta3 + 39 * zeta4 / 2 ) +
	    CA * ( - 11 * lQ23 / 9 + ( - 233. / 18. + 2 * zeta2 ) * lQ22 + ( - 2545. / 54. - 22 * zeta2 / 3 + 26 * zeta3 ) * lQ2
		   - 51157. / 648. - 337 * zeta2 / 18 + 313 * zeta3 / 9 + 22 * zeta4 ) +
	    TR * nf * ( 4 * lQ23 / 9 + 38 * lQ22 / 9 + ( 418. / 27. + 8 * zeta2 / 3 ) * lQ2
			+ 4085. / 162. + 46 * zeta2 / 9 + 4 * zeta3 / 9 ) );

      return hfct;
    };
}
