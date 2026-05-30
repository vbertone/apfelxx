//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcdqed.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/betaqed.h"
#include "apfel/betaqcdqed.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  AlphaQCDQED::AlphaQCDQED(double              const& AlphaQCDRef,
                           double              const& AlphaQEDRef,
                           double              const& MuRef,
                           std::vector<double> const& QuarkThresholds,
                           std::vector<double> const& LeptThresholds,
                           int                 const& pt,
                           int                 const& nstep):
    MatchedEvolution{matrix<double>{2, 1, {AlphaQCDRef, AlphaQEDRef}}, MuRef, ConcatenateAndSortVectors(QuarkThresholds, LeptThresholds), nstep},
                   _QuarkThresholds(QuarkThresholds),
                   _LeptThresholds(LeptThresholds),
                   _pt(pt)
  {
    // Beta function lambda function
    _BetaFunction = [=, this] (int const& nfl, matrix<double> const& a)-> matrix<double>
    {
      const std::pair<int, int> nf_nl = NFL(nfl);
      const int nf = nf_nl.first;
      const int nl = nf_nl.second;
      const matrix<double> am{2, 2, {a(0, 0), 0, 0, a(1, 0)}};
      matrix<double> bt{2, 1};

      // Leading order
      bt -= betaQCDQED(0, nf, nl) * a;

      // Next-to-leading order
      if (_pt > 0)
        bt -= am * betaQCDQED(1, nf, nl) * a;

      // Next-to-next-to-leading order (QCD only)
      if (_pt > 1)
        bt -= matrix<double>{2, 1, {pow(a(0, 0) / FourPi, 3) * beta2qcd(nf), 0}};

      // Next-to-next-to-next-to-leading order (QCD only)
      if (_pt > 2)
        bt -= matrix<double>{2, 1, {pow(a(0, 0) / FourPi, 4) * beta3qcd(nf), 0}};

      return am * bt;
    };

    // Matching condition lambda function (only acting on the QCD coupling)
    _MatchingConditions = [=, this] (bool const& Up, int const& nf, double const& Coup) -> double
    {
      const double sgn = (Up ? 1 : -1);
      const double ep  = Coup / FourPi;

      // The expression is taken from Eqs. (22) and (25) of
      // https://arxiv.org/pdf/hep-ph/0004189.pdf.
      const std::vector<double> c{
        1,
        0,
        sgn * 14. / 3.,
        sgn * pow(4, 3) *  ( 58933. / 124416. + 2. / 3. * zeta2 * ( 1. + log(2) / 3.) + 80507. / 27648. * zeta3
                             + nf * ( - 2479. / 31104. - zeta2 / 9. ) )
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
  matrix<double> AlphaQCDQED::MatchObject(bool const& Up, int const& nfl, matrix<double> const& Coup) const
  {
    const int nf = NFL(nfl).first;
    // Apply matching only if it is a quark (and not a lepton) threshold
    if ((Up ? NFL(nfl+1).first : NFL(nfl-1).first) == nf)
      return Coup;
    else
      return matrix<double> {2, 1, {_MatchingConditions(Up, (Up ? nf : nf - 1), Coup(0, 0)), Coup(1, 0)}};
  }

  //_________________________________________________________________________________
  matrix<double> AlphaQCDQED::Derivative(int const& nfl, double const&, matrix<double> const& as) const
  {
    return _BetaFunction(nfl, as);
  }

  //_________________________________________________________________________________
  matrix<double> AlphaQCDQED::betaQCDQED(int const& pt, int const& nf, int const& nl) const
  {
    matrix<double> res{2, 2};
    if (pt == 0)
      {
        res(0, 0) = beta0qcd(nf) / FourPi;
        res(1, 1) = beta0qed(nf, nl) / FourPi;
      }
    else if (pt == 1)
      {
        const double fp2 = pow(FourPi, 2);
        res(0, 0) = beta1qcd(nf) / fp2;
        res(0, 1) = beta1qcdqed(nf) / fp2;
        res(1, 0) = beta1qedqcd(nf) / fp2;
        res(1, 1) = beta1qed(nf, nl) / fp2;
      }
    else
      throw std::runtime_error(error("AlphaQCDQED::betaQCDQED", "perturbive order out of range."));

    return res;
  }

  //_________________________________________________________________________________
  std::pair<int, int> AlphaQCDQED::NFL(int const& nfl) const
  {
    const double Qr = this->_Thresholds[nfl-1] + eps8;
    return std::pair<int, int> {NF(Qr, _QuarkThresholds), NF(Qr, _LeptThresholds)};
  }
}
