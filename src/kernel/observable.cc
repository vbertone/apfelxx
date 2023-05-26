//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/observable.h"
#include "apfel/messages.h"

namespace apfel
{
  //_____________________________________________________________________________
  template<class T>
  Observable<T>::Observable(std::vector<ConvolutionPair> ConvPair):
    _ConvPair(ConvPair)
  {
    // If the convolution-pair vector is empty throw exception
    if (_ConvPair.empty())
      throw std::runtime_error(error("Observable<T>::Observable", "Vector of convolution pairs cannot be empty."));
  }

  //_____________________________________________________________________________
  template<class T>
  Observable<T>::Observable(std::function<Set<Operator>(double const&)> const& CoefficientFunctions,
                            std::function<Set<T>(double const&)>        const& Objects):
    Observable(std::vector{ConvolutionPair{CoefficientFunctions, Objects}})
  {
  }

  //_____________________________________________________________________________
  template<class T>
  void Observable<T>::AddConvolutionPair(std::function<Set<Operator>(double const&)> const& CoefficientFunctions,
                                         std::function<Set<T>(double const&)>        const& Objects)
  {
    _ConvPair.push_back(ConvolutionPair{CoefficientFunctions, Objects});
  }

  //_____________________________________________________________________________
  template<class T>
  T Observable<T>::Evaluate(double const& Q) const
  {
    Set<T> sSF = _ConvPair[0].CoefficientFunctions(Q) * _ConvPair[0].Objects(Q);
    for (int i = 1; i < (int) _ConvPair.size(); i++)
      sSF += _ConvPair[i].CoefficientFunctions(Q) * _ConvPair[i].Objects(Q);

    return sSF.Combine();
  }

  //_____________________________________________________________________________
  template<class T>
  void Observable<T>::SetObjects(std::function<Set<T>(double const&)> const& Objects, int const& ip)
  {
    if (ip < 0 && ip >= (int) _ConvPair.size())
      throw std::runtime_error(error("Observable<T>::SetObjects", "Convolution-pair index aout of range."));

    _ConvPair[ip].Objects = Objects;
  }

  //_____________________________________________________________________________
  template<class T>
  std::function<Set<Operator>(double const&)> Observable<T>::GetCoefficientFunctions(int const& ip) const
  {
    if (ip < 0 && ip >= (int) _ConvPair.size())
      throw std::runtime_error(error("Observable<T>::GetCoefficientFunctions", "Convolution-pair index aout of range."));

    return _ConvPair[ip].CoefficientFunctions;
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
