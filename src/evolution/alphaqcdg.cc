//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcdg.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/gbeta.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  AlphaQCDg::AlphaQCDg(double              const& AlphaRef,
                       double              const& MuRef,
                       std::vector<double> const& Masses,
                       std::vector<double> const& Thresholds,
                       int                 const& pt,
                       double              const& kappa):
    MatchedEvolution{AlphaRef, MuRef, Thresholds, 1},
    _pt(pt),
    _kappa(kappa)
  {
    // Matching condition lambda function
    _MatchingConditions = [=] (bool const& Up, int const& nf, double const& Coup) -> double
    {
      // Compute log of muth2 / m2
      double LogKth = 0;
      if (Masses[nf] > 0 && Thresholds[nf] > 0 && Masses[nf] != Thresholds[nf])
        LogKth = 2 * log(Thresholds[nf] / Masses[nf]);

      const double sgn = (Up ? 1 : -1);
      const double ep  = Coup / FourPi;

      // The O(as^3) matching condition does not include the
      // logarithmic terms yet. The expression is taken from Eqs. (22)
      // and (25) of https://arxiv.org/pdf/hep-ph/0004189.pdf that do
      // report the logarithmic terms instead.
      const std::vector<double> c{
        1,
        sgn * 2. / 3. * LogKth,
        4. / 9. * pow(LogKth, 2) + sgn *  38. / 3. * LogKth + sgn * 14. / 3.,
         sgn * pow(4, 3) *  ( 58933. / 124416. + 2. / 3. * zeta2 * ( 1.  + log(2) / 3.) + 80507. / 27648. * zeta3
                              + (Up ? 8941. : 8521. ) / 1728. * LogKth + (Up ? 511. : 131. ) / 576. * pow(LogKth, 2) + pow(LogKth, 3) / 216.
                              + (Up ? nf - 1 : nf) * ( - 2479. / 31104. - zeta2 / 9. - 409. / 1728. * LogKth ) )
      };
      double match = 0, powep = 1;
      for (int i = 0; i <= _pt; i++)
        {
          match += c[i] * powep;
          powep *= ep;
        }
      return Coup * match;
    };
  }

  //_________________________________________________________________________________
  AlphaQCDg::AlphaQCDg(double const& AlphaRef, double const& MuRef, std::vector<double> const& Masses, int const& pt, double const& kappa):
    AlphaQCDg{AlphaRef, MuRef, Masses, Masses, pt, kappa}
  {
  }

  //_________________________________________________________________________________
  double AlphaQCDg::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    return _MatchingConditions(Up, (Up ? nf : nf - 1), Coup);
  }

  //_________________________________________________________________________________
  double AlphaQCDg::EvolveObject(int const& nf, double const& lnmu02, double const& lnmu2, double const& as0) const
  {
    // Return immediately "as0" if "lnmu02" and "lnmu2" are equal
    if (lnmu02 == lnmu2)
      return as0;

    // Initialise evolved coupling
    const double a0 = as0 / FourPi;

    // Define lambda parameter
    const double lambda = - a0 * beta0qcd(nf) * ( 2 * log(_kappa) + lnmu2 - lnmu02 );

    // Return evolved coupling
    return as0 * ( g1beta(lambda) + a0 * ( (_pt > 0 ? g2beta(nf, _kappa, lambda) : 0) + a0 * ( (_pt > 1 ? g3beta(nf, _kappa, lambda) : 0) + a0 * (_pt > 2 ? g4beta(nf, _kappa, lambda) : 0) ) ) );
  }
}
