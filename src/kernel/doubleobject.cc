//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/doubleobject.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"
#include "apfel/tools.h"

namespace apfel
{

  //_________________________________________________________________________
  template<class T>
  DoubleObject<T>::DoubleObject(vector<term<T>> const& terms):
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
    // The product is possible only if the rhs double object contains
    // one single term.
    auto tt = o.GetTerms();
    if (tt.size() != 1)
      throw runtime_exception("DoubleObject::operator *=",
			      "The product is possible only if the r.h.s. member contains one single term (1)");

    vector<term<V>> vt;
    auto const oc  = tt[0].coefficient;
    auto const oo1 = tt[0].object1;
    auto const oo2 = tt[0].object2;
    for (auto const& t : _terms)
      vt.push_back({t.coefficient * oc, t.object1 * oo1, t.object2 * oo2});

    return DoubleObject<V>{vt};
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
  DoubleObject<T>& DoubleObject<T>::operator *= (DoubleObject<T> const& o)
  {
    // The product is possible only if the rhs double object contains
    // one single term.
    if (o.GetTerms().size() != 1)
      throw runtime_exception("DoubleObject::operator *=",
			      "The product is possible only if the r.h.s. member contains one single term (2)");

    for (auto& t : _terms)
      {
	t.coefficient *= GetTerms()[0].coefficient;
	t.object1 *= GetTerms()[0].object1;
	t.object2 *= GetTerms()[0].object2;
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
      {
	result += t.coefficient * t.object1.Evaluate(x) * t.object2.Evaluate(z);
      }
    return result;
  }

  template<>
  double DoubleObject<Operator>::Evaluate(double const&, double const&) const
  {
    throw runtime_exception("DoubleObject::Evaluate(x,z)",
			    "This function can't be used for the specialization 'Operator' of the DoubleObject class.");
  }

}
