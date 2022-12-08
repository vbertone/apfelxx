//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcdxi.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/gbeta.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  AlphaQCDxi::AlphaQCDxi(double              const& AlphaRef,
                         double              const& MuRef,
                         std::vector<double> const& Masses,
                         std::vector<double> const& Thresholds,
                         int                 const& pt,
                         double              const& xi,
                         int                 const& nstep):
    MatchedEvolution{AlphaRef, MuRef, Thresholds, nstep},
    _pt(pt),
    _xi(xi)
  {
    // Beta function lambda function.
    _BetaFunction = [=] (int const& nf, double const& as)-> double
    {
      const double lambda = - 2 * as * beta0qcd(nf) * log(_xi) / FourPi;
      double bt = 0;
      double powas = as * as;
      for (int i = 0; i <= _pt; i++)
        {
          bt += powas * betaQCD(i, nf, lambda);
          powas *= as;
        }
      return bt / 2;
    };

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
  AlphaQCDxi::AlphaQCDxi(double const& AlphaRef, double const& MuRef, std::vector<double> const& Masses, int const& pt, double const& xi, int const& nstep):
    AlphaQCDxi{AlphaRef, MuRef, Masses, Masses, pt, xi, nstep}
  {
  }

  //_________________________________________________________________________________
  double AlphaQCDxi::MatchObject(bool const& Up, int const& nf, double const& Coup) const
  {
    return _MatchingConditions(Up, nf, Coup);
  }

  //_________________________________________________________________________________
  double AlphaQCDxi::Derivative(int const& nf, double const&, double const& as) const
  {
    return _BetaFunction(nf, as);
  }

  //_________________________________________________________________________________
  double AlphaQCDxi::betaQCD(int const& pt, int const& nf, double const& lambda) const
  {
    double res;
    if (pt == 0)
      res = - 2 * beta0qcd(nf) * pow(g1beta(lambda), 2);
    else if (pt == 1)
      {
        const double bt0 = - 2 * beta0qcd(nf);
        const double b1  = - 2 * beta1qcd(nf) / bt0;
        const double g1  = g1beta(lambda);
        const double g2  = g2beta(nf, 1, lambda);
        res = bt0 * ( b1 * pow(g1, 3) + 2 * g1 * ( g2 - pow(g1, 2) * bt0 * log(_xi) ) );
      }
    else if (pt == 2)
      {
        const double bt0 = - 2 * beta0qcd(nf);
        const double b1  = - 2 * beta1qcd(nf) / bt0;
        const double b2  = - 2 * beta2qcd(nf) / bt0;
        const double g1  = g1beta(lambda);
        const double g2  = g2beta(nf, 1, lambda);
        const double g3  = g3beta(nf, 1, lambda);
        const double lx  = bt0 * log(_xi);
        res = bt0 * ( - 5 * b1 * pow(g1, 4) * lx + b2 * pow(g1, 4) + 3 * b1 * g2 * pow(g1, 2)
                      + 3 * pow(g1, 4) * pow(lx, 2) - 6 * g2 * pow(g1, 2) * lx + 2 * g3 * g1 + pow(g2, 2) );
      }
    else if (pt == 3)
      {
        const double bt0 = - 2 * beta0qcd(nf);
        const double b1  = - 2 * beta1qcd(nf) / bt0;
        const double b2  = - 2 * beta2qcd(nf) / bt0;
        const double b3  = - 2 * beta3qcd(nf) / bt0;
        const double g1  = g1beta(lambda);
        const double g2  = g2beta(nf, 1, lambda);
        const double g3  = g3beta(nf, 1, lambda);
        const double g4  = g4beta(nf, 1, lambda);
        const double lx  = bt0 * log(_xi);
        res = bt0 * ( 13 * b1 * pow(g1, 5) * pow(lx, 2) - 3 * pow(b1, 2) * pow(g1, 5) * lx - 6 * b2 * pow(g1, 5) * lx
                      - 20 * b1 * g2 * pow(g1, 3) * lx + b3 * pow(g1, 5) + 4 * b2 * g2 * pow(g1, 3) + 3 * b1 * g3 * pow(g1, 2)
                      + 3 * b1 * pow(g2, 2) * g1 - 4 * pow(g1, 5) * pow(lx, 3) + 12 * g2 * pow(g1, 3) * pow(lx, 2)
                      - 6 * g3 * pow(g1, 2) * lx - 6 * pow(g2, 2) * g1 * lx + 2 * g4 * g1 + 2 * g2 * g3 );
      }
    else
      throw std::runtime_error(error("AlphaQCDxi::betaQCD", "perturbive order out of range."));

    return res / pow(FourPi, pt + 1);
  }
}
