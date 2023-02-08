//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcd.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  AlphaQCD::AlphaQCD(double              const& AlphaRef,
                     double              const& MuRef,
                     std::vector<double> const& Masses,
                     std::vector<double> const& Thresholds,
                     int                 const& pt,
                     int                 const& nstep):
    MatchedEvolution{AlphaRef, MuRef, Thresholds, nstep},
    _pt(pt)
  {
    // Beta function lambda function
    _BetaFunction = [=] (int const& nf, double const& as)-> double
    {
      double bt = 0, powas = as * as;
      for (int i = 0; i <= _pt; i++)
        {
          bt -= powas * betaQCD(i, nf);
          powas *= as;
        }
      return bt;
    };

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
  AlphaQCD::AlphaQCD(double const& AlphaRef, double const& MuRef, std::vector<double> const& Masses, int const& pt, int const& nstep):
    AlphaQCD{AlphaRef, MuRef, Masses, Masses, pt, nstep}
  {
  }

  //_________________________________________________________________________________
  double AlphaQCD::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    return _MatchingConditions(Up, (Up ? nf : nf - 1), Coup);
  }

  //_________________________________________________________________________________
  double AlphaQCD::Derivative(int const& nf, double const&, double const& as) const
  {
    return _BetaFunction(nf, as);
  }

  //_________________________________________________________________________________
  double AlphaQCD::betaQCD(int const& pt, int const& nf) const
  {
    double res;
    if (pt == 0)
      res = beta0qcd(nf);
    else if (pt == 1)
      res = beta1qcd(nf);
    else if (pt == 2)
      res = beta2qcd(nf);
    else if (pt == 3)
      res = beta3qcd(nf);
    else
      throw std::runtime_error(error("AlphaQCD::betaQCD", "perturbive order out of range."));

    return res / pow(FourPi, pt + 1);
  }
}
