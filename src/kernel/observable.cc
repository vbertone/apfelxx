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
			 function<Set<Distribution>(double const&)> const& Distributions)
  {
    // Define observable function
    _Observable = [=] (double const& Q) -> Distribution
      {
	Set<Distribution> sSF = CoefficientFunctions(Q) * Distributions(Q);
	return sSF.Combine();
      };
  }

  //_____________________________________________________________________________
  Observable::Observable(function<Distribution(double const&)> const& Obs):
    _Observable(Obs)
  {
  }

  //_____________________________________________________________________________
  Distribution Observable::Evaluate(double const& Q) const
  {
    return _Observable(Q);
  }

}
