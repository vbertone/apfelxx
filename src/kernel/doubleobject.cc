//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/doubleobject.h"
#include "apfel/operator.h"
#include "apfel/tools.h"

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
  template<class V> DoubleObject<V> DoubleObject<T>::operator *= (DoubleObject<V> const& o) const
  {
    auto tt = o.GetTerms();
    std::vector<term<V>> vt;
    for (auto const& t : _terms)
      {
        const double tc = t.coefficient;
        for (auto const& s : tt)
          {
            const double sc = tc * s.coefficient;
            const auto o1 = t.object1 * s.object1;
            const auto o2 = t.object2 * s.object2;
            vt.push_back({sc, o1, o2});
          }
      }
    return DoubleObject<V> {vt};
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
    for (auto& t : o.GetTerms())
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
  template DoubleObject<Distribution> DoubleObject<Operator>::operator *= (DoubleObject<Distribution> const&) const;

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

  //_________________________________________________________________________
  template<>
  DoubleObject<Operator>& DoubleObject<Operator>::operator *= (DoubleObject<Operator> const& o)
  {
    auto tt = o.GetTerms();
    std::vector<term<Operator>> vt;
    for (auto const& t : _terms)
      {
        const double tc = t.coefficient;
        for (auto const& s : tt)
          {
            const double sc = tc * s.coefficient;
            const auto o1 = t.object1 * s.object1;
            const auto o2 = t.object2 * s.object2;
            vt.push_back({sc, o1, o2});
          }
      }
    // Clear "_terms" and equal it to "vt".
    _terms.clear();
    _terms = vt;
    return *this;
  }
}
