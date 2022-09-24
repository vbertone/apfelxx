//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcdg.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
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
    _lnkappa(log(kappa))
  {
    // Matching condition lambda function.
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
        sgn * pow(4, 3) *  ( - 58933. / 124416. - 2. / 3. * zeta2 * ( 1.  + log(2) / 3.) - 80507. / 27648. * zeta3 + nf * ( 2479. / 31104. + zeta2 / 9. ) )
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
    return _MatchingConditions(Up, nf, Coup);
  }

  //_________________________________________________________________________________
  double AlphaQCDg::EvolveObject(int const& nf, double const& lnmu02, double const& lnmu2, double const& as0) const
  {
    // Return immediately "as0" if "lnmu02" and "lnmu2" are equal.
    if (lnmu02 == lnmu2)
      return as0;

    // Initialise evolved coupling
    const double a0 = as0 / FourPi;

    // Define lambda parameter
    const double lambda = - a0 * beta0qcd(nf) * ( 2 * _lnkappa + lnmu2 - lnmu02 );

    // Return evolved coupling
    return as0 * ( g1beta(lambda) + a0 * ( (_pt > 0 ? g2beta(nf, lambda) : 0) + a0 * ( (_pt > 1 ? g3beta(nf, lambda) : 0) + a0 * (_pt > 2 ? g4beta(nf, lambda) : 0) ) ) );
  }

  //_________________________________________________________________________________
  double AlphaQCDg::g1beta(double const& lambda) const
  {
    return 1 / ( 1 - lambda );
  }

  //_________________________________________________________________________________
  double AlphaQCDg::g2beta(int const& nf, double const& lambda) const
  {
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1  = - 2 * beta1qcd(nf) / bt0;
    return - ( b1 * log(1 - lambda) + bt0 * _lnkappa ) / pow(1 - lambda, 2);
  }

  //_________________________________________________________________________________
  double AlphaQCDg::g3beta(int const& nf, double const& lambda) const
  {
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1  = - 2 * beta1qcd(nf) / bt0;
    const double b2  = - 2 * beta2qcd(nf) / bt0;
    const double ln1ml = log(1 - lambda);
    return - ( b2 * lambda - pow(b1, 2) * ( lambda + ln1ml - pow(ln1ml, 2) ) + bt0 * b1 * ( 2 * ln1ml - 1 ) * _lnkappa + pow(bt0 * _lnkappa, 2) ) / pow(1 - lambda, 3);
  }

  //_________________________________________________________________________________
  double AlphaQCDg::g4beta(int const& nf, double const& lambda) const
  {
    const double bt0 = - 2 * beta0qcd(nf);
    const double b1  = - 2 * beta1qcd(nf) / bt0;
    const double b2  = - 2 * beta2qcd(nf) / bt0;
    const double b3  = - 2 * beta3qcd(nf) / bt0;
    const double ln1ml = log(1 - lambda);
    return ( ( b3 - b2 * b1 ) * lambda - ( pow(b1, 3) - 2 * b2 * b1 + b3 ) * pow(lambda, 2) / 2
             + ( 2 * b1 * ( pow(b1, 2) - b2 ) * lambda - b2 * b1 ) * ln1ml
             + pow(b1, 3 ) * ( 5. / 2. - ln1ml ) * pow(ln1ml, 2)
             + ( 2 * ( pow(b1, 2) - b2 ) * lambda + pow(b1, 2) * ( 5 - 3 * ln1ml ) * ln1ml - b2 ) * bt0 * _lnkappa
             + b1 * ( - 3 * ln1ml + 5. / 2. ) * pow(bt0 * _lnkappa, 2) - pow(bt0 * _lnkappa, 3)
           ) / pow(1 - lambda, 4);
  }
}
