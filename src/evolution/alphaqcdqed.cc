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
                   _pt(pt)
  {
    // Beta function lambda function
    _BetaFunction = [=] (int const& nfl, matrix<double> const& a)-> matrix<double>
    {
      const double Qr = ConcatenateAndSortVectors(QuarkThresholds, LeptThresholds)[nfl-1] + eps8;
      const int nf = NF(Qr, QuarkThresholds);
      const int nl = NF(Qr, LeptThresholds);
      const matrix<double> am{2, 2, {a(0, 0), 0, 0, a(1, 0)}};
      matrix<double> bt{2, 1};

      // Leading order
      bt -= betaQCDQED(0, nf, nl) * a;

      // Next-to-leading order
      if (_pt > 0)
        bt -= am * betaQCDQED(1, nf, nl) * a;

      return am * bt;
    };
  }

  //_________________________________________________________________________________
  matrix<double> AlphaQCDQED::MatchObject(bool const&, int const&, matrix<double> const& Coup) const
  {
    return Coup;
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
}
