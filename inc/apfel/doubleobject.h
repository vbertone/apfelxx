//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/tools.h"

namespace apfel
{
  /**
   * @brief The term structure that contains all the objects of a
   * single term of a double object.
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
   * multiplicative constant. This mother class provides the basic
   * ingredients for the computation double convolutions required in
   * SIDIS and DY.
   */
  template<class T>
  class DoubleObject
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DoubleObject constructor.
     */
    DoubleObject();

    /**
     * @brief The DoubleObject constructor.
     * @param terms: vector of term objects of the T kind
     */
    DoubleObject(std::vector<term<T>> const& terms);
    ///@}

    /**
     * @brief Function to add more terms.
     * @param newterm: new term to be appended to the vector of terms
     */
    void AddTerm(term<T> const& newterm);

    /**
     * @brief Function to get the terms.
     * @return The vector of terms
     */
    std::vector<term<T>> GetTerms() const { return _terms; };

    /**
     * @brief Funtion that evaluates the double distribution.
     * @param x: value of the first variable
     * @param z: value of the second variable
     * @return The value of the double distribution in (x,z)
     */
    double Evaluate(double const& x, double const& z) const;

    /**
     * @name Binary operators
     */
    ///@{
    template<class V> DoubleObject<V> operator *= (DoubleObject<V> const& o) const;
    DoubleObject<T>& operator *= (double const& s);          //!< this *= scalar
    DoubleObject<T>& operator /= (double const& s);          //!< this /= scalar
    DoubleObject<T>& operator *= (DoubleObject<T> const& o); //!< this *= DoubleObject
    DoubleObject<T>& operator += (DoubleObject<T> const& o); //!< this += DoubleObject
    ///@}

  private:
    std::vector<term<T>> _terms;
  };

  /**
   * @name Ternary operators
   */
  ///@{
  template<class A, class B>
  DoubleObject<B> operator * (DoubleObject<A> lhs, DoubleObject<B> const& rhs) { return lhs *= rhs; }

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
  ///@}
}
