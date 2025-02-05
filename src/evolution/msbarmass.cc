//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/msbarmass.h"
#include "apfel/gammam.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

#include <map>

namespace apfel
{
  //_________________________________________________________________________________
  MSbarMass::MSbarMass(double                               const& MRef,
                       double                               const& MuRef,
                       std::vector<double>                  const& Masses,
                       int                                  const& pt,
                       std::function<double(double const&)> const& Alphas,
                       int                                  const& nsteps):
    MatchedEvolution{MRef, MuRef, Masses, nsteps},
    _pt(pt),
    _Alphas(Alphas)
  {
    // Collect coupling above and below threshold (only if the
    // threshold is above 1 GeV)
    std::map<int, std::pair<double, double>> AlphasTh;
    for (int nf = 1; nf <= (int) Masses.size(); nf++)
      if (Masses[nf-1] > 1)
        AlphasTh.insert({nf, std::make_pair(Alphas(Masses[nf-1] * ( 1 - eps8 )), Alphas(Masses[nf-1] * ( 1 + eps8 )))});

    // Anomalous-dimension lambda function
    _AnomalousDimension = [=] (int const& nf, double const& as)-> double
    {
      double ad = 0;
      for (int i = 0; i <= _pt; i++)
        ad -= pow(as / FourPi, i + 1) * AnomalouDimension(i, nf);
      return ad;
    };

    // Matching condition lambda function
    _MatchingConditions = [=] (bool const& Up, int const& nf, double const& mass) -> double
    {
      // Get coupling above or below threshold as appropriate
      const double cp = (Up ? AlphasTh.at(nf+1).second : AlphasTh.at(nf+1).first) / FourPi;

      // Expressions are taken from
      // https://arxiv.org/pdf/hep-ph/0004189.
      const std::vector<double> c{1, 0, - (Up ? 1 : -1) * 89. / 27.};
      double match = 0;
      for (int i = 0; i <= _pt; i++)
        match += pow(cp, i) * c[i];
      return mass * match;
    };
  }

  //_________________________________________________________________________________
  double MSbarMass::MatchObject(bool const& Up, int const& nf, double const& mass) const
  {
    return _MatchingConditions(Up, (Up ? nf : nf - 1), mass);
  }

  //_________________________________________________________________________________
  double MSbarMass::Derivative(int const& nf, double const& t, double const& mass) const
  {
    return mass * _AnomalousDimension(nf, _Alphas(exp(t / 2)));
  }

  //_________________________________________________________________________________
  double MSbarMass::AnomalouDimension(int const& pt, int const& nf) const
  {
    double res;
    if (pt == 0)
      res = gammam0();
    else if (pt == 1)
      res = gammam1(nf);
    else if (pt == 2)
      res = gammam2(nf);
    else
      throw std::runtime_error(error("MSbarMass::AnomalouDimension", "perturbive order out of range."));

    return res;
  }
}
