//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/tools.h"

#include <vector>

using std::vector;

namespace apfel
{
  /**
   * @brief The term struct
   */
  template <class V>
    struct term
    {
      double coefficient;
      V      object1;
      V      object2;
    };

  /**
   * @brief The DoubleObject class is a collection of pairs of single
   * objects (Distributions or Operators) accompained by a
   * multiplicative constant.
   *
   * This mother class provides the basic ingredients for the
   * computation double convolutions required in SIDIS and DY.
   */
  template<class T>
  class DoubleObject
  {
  public:
    /**
     * @brief DoubleObject constructor
     */
    DoubleObject(vector<term<T>> const& terms);

    /**
     * @brief AddTerm
     */
    void AddTerm(term<T> const& newterm);

    /**
     * @brief GetTerms
     */
    vector<term<T>> GetTerms() const { return _terms; };

    /**
     * @brief Evaluate
     */
    double Evaluate(double const& x, double const& z) const;

    /**
     * @brief operators *= and /=
     */
    template<class V> DoubleObject<V> operator *= (DoubleObject<V> const& o) const;
    DoubleObject<T>& operator *= (double const& s);          //!< this *= scalar
    DoubleObject<T>& operator /= (double const& s);          //!< this /= scalar
    DoubleObject<T>& operator *= (DoubleObject<T> const& o); //!< this *= DoubleObject
    DoubleObject<T>& operator += (DoubleObject<T> const& o); //!< this += DoubleObject

  private:
    vector<term<T>> _terms;
  };

  /**
   * @brief operator * and / definition
   */
  template<class A, class B>
  DoubleObject<B> operator * (DoubleObject<A> lhs, DoubleObject<B> const& rhs) { return lhs *= rhs; }

  // other operators
  template<class T>
  DoubleObject<T> operator * (double const& s, DoubleObject<T> rhs) { return rhs *= s; }

  template<class T>
  DoubleObject<T> operator * (DoubleObject<T> lhs, double const& s) { return lhs *= s; }

  template<class T>
  DoubleObject<T> operator / (double const& s, DoubleObject<T> rhs) { return rhs /= s; }

  template<class T>
  DoubleObject<T> operator / (DoubleObject<T> lhs, double const& s) { return lhs /= s; }

  template<class T>
  DoubleObject<T> operator * (DoubleObject<T> lhs, DoubleObject<T> const& rhs) { return lhs *= rhs; }

  template<class T>
  DoubleObject<T> operator + (DoubleObject<T> lhs, DoubleObject<T> const& rhs) { return lhs += rhs; }

}
