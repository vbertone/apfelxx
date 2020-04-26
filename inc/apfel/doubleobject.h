//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/tools.h"

#include <functional>

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
     * @brief Function that evaluates the double distribution.
     * @param x: value of the first variable
     * @param z: value of the second variable
     * @return The value of the double distribution in (x, z)
     */
    double Evaluate(double const& x, double const& z) const;

    /**
     * @brief Function that evaluates the double distribution in one of
     * the two variables leaving the other undetermined.
     * @param iv: variable index. It can be either 1 or 2.
     * @param y: value of the iv-th variable
     * @return The value of the double object in y for any other value
     * of the other variable
     */
    T Evaluate(int const& iv, double const& y) const;

    /**
     * @brief Function that evaluates the derivative of the double
     * distribution.
     * @param x: value of the first variable
     * @param z: value of the second variable
     * @return The value of the derivative of the double distribution
     * in (x, z)
     */
    double Derive(double const& x, double const& z) const;

    /**
     * @brief Function that evaluates the derivative of a double
     * distribution in one of the two variables leaving the other
     * undetermined.
     * @param iv: variable index. It can be either 1 or 2.
     * @param y: value of the iv-th variable in which the derivative
     * has to be computed
     * @return The value of the derivative of the double object in y
     * for any other value of the other variable
     */
    T Derive(int const& iv, double const& y) const;

    /**
     * @brief Function that evaluates the integral of the double
     * distribution.
     * @param xl: value of the lower bound of the of the first variable
     * @param xu: value of the upper bound of the of the first variable
     * @param zl: value of the lower bound of the of the second variable
     * @param zu: value of the upper bound of the of the second variable
     * @return The value of the integral of the double distribution
     */
    double Integrate(double const& xl, double const& xu, double const& zl, double const& zu) const;

    /**
     * @brief Function that evaluates the integral of the double
     * distribution in one of the two variables leaving the other
     * undetermined.
     * @param iv: variable index. It can be either 1 or 2.
     * @param yl: value of the lower bound of the of the iv-th variable
     * @param yu: value of the upper bound of the of the iv-th variable
     * @return The value of the integral of the double distribution in
     * [yl: yu] for any other value of the other variable
     */
    T Integrate(int const& iv, double const& yl, double const& yu) const;

    /**
     * @brief This function multiplies the distributions single terms
     * of the DoubleObject by a respective function.
     * @param fx: that function that multiplies the first distribution
     * @param fz: that function that multiplies the second distribution
     */
    DoubleObject<T>& MultiplyBy(std::function<double(double const&)> const& fx, std::function<double(double const&)> const& fz);

    /**
     * @name Binary operators
     */
    ///@{
    template<class V> DoubleObject<V> operator *= (DoubleObject<V> const& o) const;
    DoubleObject<T>& operator *= (double const& s);          //!< this *= scalar
    DoubleObject<T>& operator *= (DoubleObject<T> const& o); //!< this *= DoubleObject
    DoubleObject<T>& operator /= (double const& s);          //!< this /= scalar
    DoubleObject<T>& operator += (DoubleObject<T> const& o); //!< this += DoubleObject
    DoubleObject<T>& operator -= (DoubleObject<T> const& o); //!< this -= DoubleObject
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

  template<class T>
  DoubleObject<T> operator - (DoubleObject<T> lhs, DoubleObject<T> const& rhs) { return lhs -= rhs; }
  ///@}
}
