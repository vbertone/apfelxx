//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
  template<class T, class U>
  DoubleObject<T, U>::DoubleObject():
    _terms({})
  {
  }

  //_________________________________________________________________________
  template<class T, class U>
  DoubleObject<T, U>::DoubleObject(std::vector<term<T, U>> const& terms):
    _terms(terms)
  {
  }

  //_________________________________________________________________________
  template<class T, class U>
  void DoubleObject<T, U>::AddTerm(term<T, U> const& newterm)
  {
    _terms.push_back(newterm);
  }

  //_________________________________________________________________________
  template<class T, class U>
  template<class V> DoubleObject<V> DoubleObject<T, U>::operator *= (DoubleObject<V> const& o) const
  {
    const std::vector<term<V>> tt = o.GetTerms();
    std::vector<term<V>> vt;
    for (auto const& t : _terms)
      {
        const double tc = t.coefficient;
        for (auto const& s : tt)
          {
            const double sc = tc * s.coefficient;
            const V o1 = t.object1 * s.object1;
            const V o2 = t.object2 * s.object2;
            vt.push_back({sc, o1, o2});
          }
      }
    return DoubleObject<V> {vt};
  }

  //_________________________________________________________________________
  template<class T, class U>
  DoubleObject<T, U>& DoubleObject<T, U>::operator *= (std::function<double(double const&)> const& f)
  {
    for (auto& t : _terms)
      {
        t.object1 *= f;
        t.object2 *= f;
      }
    return *this;
  }

  //_________________________________________________________________________
  template<class T, class U>
  DoubleObject<T, U>& DoubleObject<T, U>::operator *= (double const& s)
  {
    for (auto& t : _terms)
      t.coefficient *= s;

    return *this;
  }

  //_________________________________________________________________________
  template<class T, class U>
  DoubleObject<T, U>& DoubleObject<T, U>::operator /= (double const& s)
  {
    for (auto& t : _terms)
      t.coefficient /= s;

    return *this;
  }

  //_________________________________________________________________________
  template<class T, class U>
  DoubleObject<T, U>& DoubleObject<T, U>::operator += (DoubleObject<T, U> const& o)
  {
    for (auto const& t : o.GetTerms())
      _terms.push_back(t);

    return *this;
  }

  //_________________________________________________________________________
  template<class T, class U>
  DoubleObject<T, U>& DoubleObject<T, U>::operator -= (DoubleObject<T, U> const& o)
  {
    for (auto& t : o.GetTerms())
      {
        t.coefficient *= -1;
        _terms.push_back(t);
      }
    return *this;
  }

  //_________________________________________________________________________
  template<class V, class W>
  std::ostream& operator << (std::ostream& os, DoubleObject<V, W> const& dob)
  {
    const std::vector<term<V, W>> tms = dob.GetTerms();
    os << "DoubleObject: " << &dob << "\n";
    for (int i = 0; i < (int) tms.size(); i++)
      {
        os << "Term: " << i << "\n";
        os << "- Coeffient: " << tms[i].coefficient << "\n";
        os << "- Object1:\n" << tms[i].object1 << "\n";
        os << "- Object2:\n" << tms[i].object2 << "\n";
      }
    return os;
  }

  // Specialisations
  //_________________________________________________________________________________
  template class DoubleObject<Distribution>;
  template class DoubleObject<Operator>;
  template class DoubleObject<Operator, Distribution>;
  template class DoubleObject<Distribution, Operator>;
  template DoubleObject<Distribution> DoubleObject<Operator>::operator *= (DoubleObject<Distribution> const&) const;
  template std::ostream& operator << (std::ostream& os, DoubleObject<Distribution> const& dob);
  template std::ostream& operator << (std::ostream& os, DoubleObject<Operator> const& dob);
  template std::ostream& operator << (std::ostream& os, DoubleObject<Distribution, Operator> const& dob);
  template std::ostream& operator << (std::ostream& os, DoubleObject<Operator, Distribution> const& dob);

  //_________________________________________________________________________________
  template<>
  DoubleObject<Distribution>& DoubleObject<Distribution>::MultiplyBy(std::function<double(double const&)> const& fx, std::function<double(double const&)> const& fz)
  {
    // Return immediately if "_terms" is empty
    if (_terms.empty())
      return *this;

    // Take the grids from the first element of the vector of terms
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
      result += t.coefficient * t.object1.Evaluate(x) * t.object2.Evaluate(z);

    return result;
  }

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Evaluate1(double const& x) const
  {
    Distribution result = _terms[0].coefficient * _terms[0].object1.Evaluate(x) * _terms[0].object2;
    for (int i = 1; i < (int) _terms.size(); i++)
      result += _terms[i].coefficient * _terms[i].object1.Evaluate(x) * _terms[i].object2;

    return result;
  }

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Evaluate2(double const& z) const
  {
    Distribution result = _terms[0].coefficient * _terms[0].object2.Evaluate(z) * _terms[0].object1;
    for (int i = 1; i < (int) _terms.size(); i++)
      result += _terms[i].coefficient * _terms[i].object2.Evaluate(z) * _terms[i].object1;

    return result;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Derive(double const& x, double const& z) const
  {
    double result = 0;
    for (auto const& t : _terms)
      result += t.coefficient * t.object1.Derive(x) * t.object2.Derive(z);

    return result;
  }

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Derive1(double const& x) const
  {
    Distribution result = _terms[0].coefficient * _terms[0].object1.Derive(x) * _terms[0].object2;
    for (int i = 1; i < (int) _terms.size(); i++)
      result += _terms[i].coefficient * _terms[i].object1.Derive(x) * _terms[i].object2;

    return result;
  }

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Derive2(double const& z) const
  {
    Distribution result = _terms[0].coefficient * _terms[0].object2.Derive(z) * _terms[0].object1;
    for (int i = 1; i < (int) _terms.size(); i++)
      result += _terms[i].coefficient * _terms[i].object2.Derive(z) * _terms[i].object1;

    return result;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Integrate(double const& xl, double const& xu, double const& zl, double const& zu) const
  {
    double result = 0;
    for (auto const& t : _terms)
      result += t.coefficient * t.object1.Integrate(xl, xu) * t.object2.Integrate(zl, zu);

    return result;
  }

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Integrate1(double const& xl, double const& xu) const
  {
    Distribution result = _terms[0].coefficient * _terms[0].object1.Integrate(xl, xu) * _terms[0].object2;
    for (int i = 1; i < (int) _terms.size(); i++)
      result += _terms[i].coefficient * _terms[i].object1.Integrate(xl, xu) * _terms[i].object2;

    return result;
  }

  //_________________________________________________________________________________
  template<>
  Distribution DoubleObject<Distribution>::Integrate2(double const& zl, double const& zu) const
  {
    Distribution result = _terms[0].coefficient * _terms[0].object2.Integrate(zl, zu) * _terms[0].object1;
    for (int i = 1; i < (int) _terms.size(); i++)
      result += _terms[i].coefficient * _terms[i].object2.Integrate(zl, zu) * _terms[i].object1;

    return result;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Integrate(double const& xl, double const& xu, std::function<double(double const&)> zlx, std::function<double(double const&)> zux) const
  {
    double result = 0;
    for (auto const& t : _terms)
      // Multiply the first object by the function of the first
      // variable given by the integral of the second object over the
      // x-dependent range. Finally integrate over x.
      result += t.coefficient * ([=] (double const& x) -> double{ return t.object2.Integrate(zlx(x), zux(x)); } * t.object1).Integrate(xl, xu);

    return result;
  }

  //_________________________________________________________________________________
  template<>
  double DoubleObject<Distribution>::Integrate(std::function<double(double const&)> xlz, std::function<double(double const&)> xuz, double const& zl, double const& zu) const
  {
    double result = 0;
    for (auto const& t : _terms)
      // Multiply the second object by the function of the second
      // variable given by the integral of the first object over the
      // z-dependent range. Finally integrate over z.
      result += t.coefficient * ([=] (double const& z) -> double{ return t.object1.Integrate(xlz(z), xuz(z)); } * t.object2).Integrate(zl, zu);

    return result;
  }

  //_________________________________________________________________________
  template<>
  DoubleObject<Operator>& DoubleObject<Operator>::operator *= (DoubleObject<Operator> const& o)
  {
    std::vector<term<Operator>> tt = o.GetTerms();
    std::vector<term<Operator>> vt;
    for (auto const& t : _terms)
      {
        const double tc = t.coefficient;
        for (auto const& s : tt)
          {
            const double sc = tc * s.coefficient;
            const Operator o1 = t.object1 * s.object1;
            const Operator o2 = t.object2 * s.object2;
            vt.push_back({sc, o1, o2});
          }
      }
    // Clear "_terms" and equal it to "vt".
    _terms.clear();
    _terms = vt;

    return *this;
  }
}
