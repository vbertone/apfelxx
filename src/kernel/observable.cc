//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/observable.h"

namespace apfel {
  //_____________________________________________________________________________
  Observable::Observable(function<Set<Operator>(double const&)>     const& CoefficientFunctions,
			 function<Set<Distribution>(double const&)> const& Distributions):
    _CoefficientFunctions(CoefficientFunctions),
    _Distributions(Distributions)
  {
  }

  //_____________________________________________________________________________
  Distribution Observable::Evaluate(double const& Q) const
  {
    const auto sSF = _CoefficientFunctions(Q) * _Distributions(Q);
    return sSF.Combine();
  }

  //_____________________________________________________________________________
  double Observable::Evaluate(double const& x, double const& Q) const
  {
    return this->Evaluate(Q).Evaluate(x);
  }
}
