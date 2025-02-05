//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcdmsbarmass.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  AlphaQCDMSbarMass::AlphaQCDMSbarMass(double              const& AlphaRef,
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

      // The expression is taken from Eqs. (2.20) of
      // https://arxiv.org/pdf/1605.01946.
      const std::vector<double> c{
        1,
        sgn * 2. / 3. * LogKth,
        4. / 9. * pow(LogKth, 2) + sgn *  22. / 3. * LogKth - sgn * 22. / 9.
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
  AlphaQCDMSbarMass::AlphaQCDMSbarMass(double const& AlphaRef, double const& MuRef, std::vector<double> const& Masses, int const& pt, int const& nstep):
    AlphaQCDMSbarMass{AlphaRef, MuRef, Masses, Masses, pt, nstep}
  {
  }

  //_________________________________________________________________________________
  double AlphaQCDMSbarMass::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    return _MatchingConditions(Up, (Up ? nf : nf - 1), Coup);
  }

  //_________________________________________________________________________________
  double AlphaQCDMSbarMass::Derivative(int const& nf, double const&, double const& as) const
  {
    return _BetaFunction(nf, as);
  }

  //_________________________________________________________________________________
  double AlphaQCDMSbarMass::betaQCD(int const& pt, int const& nf) const
  {
    double res;
    if (pt == 0)
      res = beta0qcd(nf);
    else if (pt == 1)
      res = beta1qcd(nf);
    else if (pt == 2)
      res = beta2qcd(nf);
    else
      throw std::runtime_error(error("AlphaQCDMSbarMass::betaQCD", "perturbive order out of range."));

    return res / pow(FourPi, pt + 1);
  }
}
