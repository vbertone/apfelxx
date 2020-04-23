//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubleobject.h"
#include "apfel/operator.h"
#include "apfel/tools.h"
#include "apfel/messages.h"
#include "apfel/constants.h"

#include <cmath>

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>::DoubleObject():
    _terms({})
  {
  }

  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>::DoubleObject(std::vector<term<T>> const& terms):
    _terms(terms)
  {
  }

  //_________________________________________________________________________
  template<class T>
  void DoubleObject<T>::AddTerm(term<T> const& newterm)
  {
    _terms.push_back(newterm);
  }

  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>& DoubleObject<T>::operator *= (double const& s)
  {
    for (auto& t : _terms)
      t.coefficient *= s;

    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>& DoubleObject<T>::operator /= (double const& s)
  {
    for (auto& t : _terms)
      t.coefficient /= s;

    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>& DoubleObject<T>::operator += (DoubleObject<T> const& o)
  {
    for (auto const& t : o.GetTerms())
      _terms.push_back(t);

    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>& DoubleObject<T>::operator -= (DoubleObject<T> const& o)
  {
    for (auto& t : o.GetTerms())
      {
        t.coefficient *= -1;
        _terms.push_back(t);
      }

    return *this;
  }

  // Specializations
  //_________________________________________________________________________________
  template class DoubleObject<Distribution>;
  template class DoubleObject<Operator>;

  //_________________________________________________________________________________
  template<>
  DoubleObject<Distribution>& DoubleObject<Distribution>::MultiplyBy(std::function<double(double const&)> const& fx, std::function<double(double const&)> const& fz)
  {
    // Take the grids from the first element of the vector of terms.
    const Distribution dfx{_terms[0].object1.GetGrid(), fx};
    const Distribution dfz{_terms[0].object2.GetGrid(), fz};
    for (auto& t : _terms)
      {
        t.object1 *= dfx;
        t.object2 *= dfz;
      }

    return *this;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Evaluate(double const& x, double const& z) const
  {
    double result = 0;
    for (auto const& t : _terms)
      if (t.coefficient == 1)
        result += t.object1.Evaluate(x) * t.object2.Evaluate(z);
      else
        result += t.coefficient * t.object1.Evaluate(x) * t.object2.Evaluate(z);
    return result;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Derive(double const& x, double const& z) const
  {
    double result = 0;
    for (auto const& t : _terms)
      if (t.coefficient == 1)
        result += t.object1.Derive(x) * t.object2.Derive(z);
      else
        result += t.coefficient * t.object1.Derive(x) * t.object2.Derive(z);
    return result;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Integrate(double const& xl, double const& xu, double const& zl, double const& zu) const
  {
    double result = 0;
    for (auto const& t : _terms)
      if (t.coefficient == 1)
        result += t.object1.Integrate(xl, xu) * t.object2.Integrate(zl, zu);
      else
        result += t.coefficient * t.object1.Integrate(xl, xu) * t.object2.Integrate(zl, zu);
    return result;
  }
}
