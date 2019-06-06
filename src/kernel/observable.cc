//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/observable.h"

namespace apfel
{
  //_____________________________________________________________________________
  template<class T>
  Observable<T>::Observable(std::function<Set<Operator>(double const&)>     const& CoefficientFunctions,
                            std::function<Set<T>(double const&)> const& Objects):
    _CoefficientFunctions(CoefficientFunctions),
    _Objects(Objects)
  {
  }

  //_____________________________________________________________________________
  template<class T>
  T Observable<T>::Evaluate(double const& Q) const
  {
    const Set<T> sSF = _CoefficientFunctions(Q) * _Objects(Q);
    return sSF.Combine();
  }

  // Specializations
  //_________________________________________________________________________________
  template class Observable<Distribution>;
  template class Observable<Operator>;

  //_____________________________________________________________________________
  template<>
  double Observable<Distribution>::Evaluate(double const& x, double const& Q) const
  {
    return this->Evaluate(Q).Evaluate(x);
  }
}
