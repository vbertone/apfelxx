//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/alphaqcdlambda.h"
#include "apfel/constants.h"
#include "apfel/betaqcd.h"
#include "apfel/tools.h"

#include <iostream>

namespace apfel
{
  //_________________________________________________________________________________
  AlphaQCDLambda::AlphaQCDLambda(double              const& LambdaQCD,
                                 int                 const& nfRef,
                                 std::vector<double> const& Thresholds,
                                 int                 const& pt):
    _Thresholds(Thresholds),
    _pt(pt)
  {
    // Function to compute LambdaQCD at nl = nf - 1 give LambdaQCD at
    // nf. Given Eq. (11) of https://arxiv.org/pdf/hep-ph/9706430.
    const std::function<double(double const&, int const&)> LambdaQCDnl = [=] (double const& LambdaQCDnf, int const& nf) -> double
    {
      // Useful definitions
      const double nl   = nf - 1;
      const double Lh   = 2 * log(_Thresholds[nf] / _LambdaQCD[nf]);
      const double lnLh = log(Lh);

      // Relevant coefficients
      const double beta0p = beta0qcd(nl);
      const double b1p    = beta1qcd(nl) / beta0p;
      const double b2p    = beta2qcd(nl) / beta0p;
      const double b3p    = beta3qcd(nl) / beta0p;

      const double beta0 = beta0qcd(nf);
      const double b1    = beta1qcd(nf) / beta0;
      const double b2    = beta2qcd(nf) / beta0;
      const double b3    = beta3qcd(nf) / beta0;

      // Coeffients C2 and C3 given in Eq. (10) of
      // https://arxiv.org/pdf/hep-ph/9706430.
      const double C2 = - 14. / 3.;
      const double C3 = pow(4, 3) * ( - 80507. / 27648. * zeta3 - 2. / 3. * zeta2 * ( 1. + log(2) / 3.) - 58933. / 124416. + nl * ( zeta2 + 2479. / 3456. ) / 9. );

      // Leading order
      double LambdaQCDnl = LambdaQCDnf * exp( ( beta0p - beta0 ) * Lh / 2. / beta0p );

      // Next-to-leading order
      if (_pt > 0)
        LambdaQCDnl *= exp( ( ( b1p - b1 ) * lnLh - b1p * log(beta0p / beta0) ) / 2. / beta0p );

      // Next-to-next-to-leading order
      if (_pt > 1)
        LambdaQCDnl *= exp( ( ( b1p - b1 ) * b1 * lnLh + pow(b1p, 2) - pow(b1, 2) - b2p + b2 + C2 ) / ( beta0 * Lh ) / 2. / beta0p );

      // Next-to-next-to-next-to-leading order
      if (_pt > 2)
        LambdaQCDnl *= exp( ( - ( b1p - b1 ) * pow(b1 * lnLh, 2) / 2 + b1 * ( - b1p * ( b1p - b1 ) + b2p - b2 - C2 ) * lnLh
                              + ( - pow(b1p, 3) - pow(b1, 3) - b3p + b3 ) / 2 + b1p * ( pow(b1, 2) + b2p - b2 - C2 ) + C3 ) / pow(beta0 * Lh, 2) / 2. / beta0p );

      return LambdaQCDnl;
    };

    // Function to compute LambdaQCD at nf = nl + 1 give LambdaQCD at
    // nl. Derived from Eq. (11) of https://arxiv.org/pdf/hep-ph/9706430.
    const std::function<double(double const&, int const&)> LambdaQCDnf = [=] (double const& LambdaQCDnl, int const& nl) -> double
    {
      // Useful definitions
      const double nf   = nl + 1;
      const double Lh   = 2 * log(_Thresholds[nf] / _LambdaQCD[nl]);
      const double lnLh = log(Lh);

      // Relevant coefficients
      const double beta0p = beta0qcd(nf);
      const double b1p    = beta1qcd(nf) / beta0p;
      const double b2p    = beta2qcd(nf) / beta0p;
      const double b3p    = beta3qcd(nf) / beta0p;

      const double beta0 = beta0qcd(nl);
      const double b1    = beta1qcd(nl) / beta0;
      const double b2    = beta2qcd(nl) / beta0;
      const double b3    = beta3qcd(nl) / beta0;

      // Coeffients C2 and C3 given in Eq. (10) of
      // https://arxiv.org/pdf/hep-ph/9706430.
      const double C2 = 14. / 3.;
      const double C3 = - pow(4, 3) * ( - 80507. / 27648. * zeta3 - 2. / 3. * zeta2 * ( 1. + log(2) / 3.) - 58933. / 124416. + nl * ( zeta2 + 2479. / 3456. ) / 9. );

      // Leading order
      double LambdaQCDnf = LambdaQCDnl * exp( ( beta0p - beta0 ) * Lh / 2. / beta0p );

      // Next-to-leading order
      if (_pt > 0)
        LambdaQCDnf *= exp( ( ( b1p - b1 ) * lnLh - b1p * log(beta0p / beta0) ) / 2. / beta0p );

      // Next-to-next-to-leading order
      if (_pt > 1)
        LambdaQCDnf *= exp( ( ( b1p - b1 ) * b1 * lnLh + pow(b1p, 2) - pow(b1, 2) - b2p + b2 + C2 ) / ( beta0 * Lh ) / 2. / beta0p );

      // Next-to-next-to-next-to-leading order
      if (_pt > 2)
        LambdaQCDnf *= exp( ( - ( b1p - b1 ) * pow(b1 * lnLh, 2) / 2 + b1 * ( - b1p * ( b1p - b1 ) + b2p - b2 - C2 ) * lnLh
                              + ( - pow(b1p, 3) - pow(b1, 3) - b3p + b3 ) / 2 + b1p * ( pow(b1, 2) + b2p - b2 - C2 ) + C3 ) / pow(beta0 * Lh, 2) / 2. / beta0p );

      return LambdaQCDnf;
    };

    // Resize vector of LambdaQCD according to Thresholds
    _LambdaQCD.resize(_Thresholds.size());

    // Place reference value
    _LambdaQCD[nfRef - 1] = LambdaQCD;

    // Now compute LambdaQCD for nf < nfRef ...
    for (int i = nfRef - 2; i >= 0; i--)
      if (_Thresholds[i + 1] <= 0)
        _LambdaQCD[i] = 0;
      else
        _LambdaQCD[i] = LambdaQCDnl(_LambdaQCD[i + 1], i + 1);

    // ... and for nf > nfRef
    for (int i = nfRef; i < (int) _Thresholds.size(); i++)
      if (_Thresholds[i] <= 0)
        _LambdaQCD[i] = 0;
      else
        _LambdaQCD[i] = LambdaQCDnf(_LambdaQCD[i - 1], i - 1);
  }

  //_________________________________________________________________________________
  std::complex<double> AlphaQCDLambda::Evaluate(std::complex<double> const& mu2) const
  {
    // Implementation of Eq. (3) of
    // https://arxiv.org/pdf/hep-ph/9706430

    // Get number of active flavours
    const int nf = NF(sqrt(sqrt(norm(mu2))), _Thresholds);

    // Get log of LambdaQCD with the correct number of active flavours
    const std::complex<double> L = log(mu2 / pow(_LambdaQCD[nf - 1], 2));

    // Get value of beta0 at nf
    const double beta0 = beta0qcd(nf);

    // Leading order
    std::complex<double> a = 1. / ( L * beta0 );

    // Next-to-leading order
    if (_pt > 0)
      {
        const std::complex<double> lnL = log(L);
        const double b1 = beta1qcd(nf) / beta0;
        a += - b1 * lnL / pow(L * beta0, 2);

        // Next-to-next-to-leading order
        if (_pt > 1)
          {
            const double b2 = beta2qcd(nf) / beta0;
            a += ( pow(b1, 2) * ( pow(lnL, 2) - lnL - 1. ) + b2 ) / pow(L * beta0, 3);

            // Next-to-next-to-next-to-leading order
            if (_pt > 2)
              {
                const double b3 = beta3qcd(nf) / beta0;
                a += ( pow(b1, 3) * ( - pow(lnL, 3) + 5. * pow(lnL, 2) / 2. + 2. * lnL - 1. / 2. ) - 3 * b1 * b2 * lnL + b3 / 2 ) / pow(L * beta0, 4);
              }
          }
      }

    return FourPi * a;
  }

  //_________________________________________________________________________________
  double AlphaQCDLambda::Evaluate(double const& mu) const
  {
    return Evaluate(std::complex<double> {pow(mu, 2), 0}).real();
  }
}
