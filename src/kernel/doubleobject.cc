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
  template DoubleObject<Distribution> DoubleObject<Operator>::operator *= (DoubleObject<Distribution> const&) const;

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
  Distribution DoubleObject<Distribution>::Evaluate(int const& iv, double const& y) const
  {
    if (iv < 1 || iv >2)
      throw std::runtime_error(error("Evaluate", "Function index out of range: it can be either 1 or 2."));

    if (iv == 1)
      {
        Distribution result = _terms[0].coefficient * _terms[0].object1.Evaluate(y) * _terms[0].object2;
        for (int i = 1; i < (int) _terms.size(); i++)
          result += _terms[i].coefficient * _terms[i].object1.Evaluate(y) * _terms[i].object2;
        return result;
      }
    else
      {
        Distribution result = _terms[0].coefficient * _terms[0].object2.Evaluate(y) * _terms[0].object1;
        for (int i = 1; i < (int) _terms.size(); i++)
          result += _terms[i].coefficient * _terms[i].object2.Evaluate(y) * _terms[i].object1;
        return result;
      }
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
  Distribution DoubleObject<Distribution>::Derive(int const& iv, double const& y) const
  {
    if (iv < 1 || iv >2)
      throw std::runtime_error(error("Derive", "Function index out of range: it can be either 1 or 2."));

    if (iv == 1)
      {
        Distribution result = _terms[0].coefficient * _terms[0].object1.Derive(y) * _terms[0].object2;
        for (int i = 1; i < (int) _terms.size(); i++)
          result += _terms[i].coefficient * _terms[i].object1.Derive(y) * _terms[i].object2;
        return result;
      }
    else
      {
        Distribution result = _terms[0].coefficient * _terms[0].object2.Derive(y) * _terms[0].object1;
        for (int i = 1; i < (int) _terms.size(); i++)
          result += _terms[i].coefficient * _terms[i].object2.Derive(y) * _terms[i].object1;
        return result;
      }
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

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Integrate(int const& iv, double const& yl, double const& yu) const
  {
    if (iv < 1 || iv >2)
      throw std::runtime_error(error("Integrate", "Function index out of range: it can be either 1 or 2."));

    if (iv == 1)
      {
        Distribution result = _terms[0].coefficient * _terms[0].object1.Integrate(yl, yu) * _terms[0].object2;
        for (int i = 1; i < (int) _terms.size(); i++)
          result += _terms[i].coefficient * _terms[i].object1.Integrate(yl, yu) * _terms[i].object2;
        return result;
      }
    else
      {
        Distribution result = _terms[0].coefficient * _terms[0].object2.Integrate(yl, yu) * _terms[0].object1;
        for (int i = 1; i < (int) _terms.size(); i++)
          result += _terms[i].coefficient * _terms[i].object2.Integrate(yl, yu) * _terms[i].object1;
        return result;
      }
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
