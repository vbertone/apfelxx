//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
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
  template <class V, class U = V>
  struct term
  {
    double coefficient;
    V      object1;
    U      object2;
  };

  /**
   * @brief The DoubleObject class is a collection of pairs of single
   * objects (Distributions or Operators) accompained by a
   * multiplicative constant. This mother class provides the basic
   * ingredients for the computation double convolutions required in
   * SIDIS and DY.
   */
  template<class T, class U = T>
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
    DoubleObject(std::vector<term<T, U>> const& terms);
    ///@}

    /**
     * @brief Function to add more terms.
     * @param newterm: new term to be appended to the vector of terms
     */
    void AddTerm(term<T, U> const& newterm);

    /**
     * @brief Function to get the terms.
     * @return The vector of terms
     */
    std::vector<term<T, U>> GetTerms() const { return _terms; };

    /**
     * @brief Function that evaluates the double distribution.
     * @param x: value of the first variable
     * @param z: value of the second variable
     * @return The value of the double distribution in (x, z)
     */
    double Evaluate(double const& x, double const& z) const;

    /**
     * @brief Function that evaluates the double object in the first
     * variable leaving the second undetermined.
     * @param x: value of the first variable
     */
    T Evaluate1(double const& x) const;

    /**
     * @brief Function that evaluates the double object in the second
     * variable leaving the first undetermined.
     * @param z: value of the second variable
     */
    U Evaluate2(double const& z) const;

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
     * @brief Function that evaluates the derivative of the double
     * object in the first variable leaving the second undetermined.
     * @param x: value of the first variable
     */
    T Derive1(double const& x) const;

    /**
     * @brief Function that evaluates the derivative of the double
     * object in the second variable leaving the first undetermined.
     * @param z: value of the second variable
     */
    U Derive2(double const& z) const;

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
     * object in the first variable leaving the second undetermined.
     * @param xl: value of the lower bound of the of the first variable
     * @param xu: value of the upper bound of the of the first variable
     */
    T Integrate1(double const& xl, double const& xu) const;

    /**
     * @brief Function that evaluates the derivative of the double
     * object in the second variable leaving the first undetermined.
     * @param zl: value of the lower bound of the of the second variable
     * @param zu: value of the upper bound of the of the second variable
     */
    U Integrate2(double const& zl, double const& zu) const;

    /**
     * @brief Function that evaluates the integral of the double
     * distribution.
     * @param xl: value of the lower bound of the of the first variable
     * @param xu: value of the upper bound of the of the first variable
     * @param zlx: function that delimits the lower bound of the integral in z as a function of x
     * @param zux: function that delimits the upper bound of the integral in z as a function of x
     * @return The value of the integral of the double distribution
     */
    double Integrate(double const& xl, double const& xu, std::function<double(double const&)> zlx, std::function<double(double const&)> zux) const;

    /**
     * @brief Function that evaluates the integral of the double
     * distribution.
     * @param xlz: function that delimits the lower bound of the integral in x as a function of z
     * @param xuz: function that delimits the upper bound of the integral in x as a function of z
     * @param zl: value of the lower bound of the of the second variable
     * @param zu: value of the upper bound of the of the seconf variable
     * @return The value of the integral of the double distribution
     */
    double Integrate(std::function<double(double const&)> xlz, std::function<double(double const&)> xuz, double const& zl, double const& zu) const;

    /**
     * @brief This function multiplies the objects of the single terms
     * of the DoubleObject by a respective function.
     * @param fx: that function that multiplies the first distribution
     * @param fz: that function that multiplies the second distribution
     */
    DoubleObject<T, U>& MultiplyBy(std::function<double(double const&)> const& fx, std::function<double(double const&)> const& fz);

    /**
     * @name Binary operators
     */
    ///@{
    template<class V> DoubleObject<V> operator *= (DoubleObject<V> const& o) const;
    DoubleObject<T, U>& operator *= (double const& s);                               //!< this *= scalar
    DoubleObject<T, U>& operator *= (DoubleObject<T, U> const& o);                   //!< this *= DoubleObject
    DoubleObject<T, U>& operator *= (std::function<double(double const&)> const& f); //!< this *= Function of the integration variable
    DoubleObject<T, U>& operator /= (double const& s);                               //!< this /= scalar
    DoubleObject<T, U>& operator += (DoubleObject<T, U> const& o);                   //!< this += DoubleObject
    DoubleObject<T, U>& operator -= (DoubleObject<T, U> const& o);                   //!< this -= DoubleObject
    ///@}

  private:
    std::vector<term<T, U>> _terms;

    template<class V, class W>
    friend std::ostream& operator << (std::ostream& os, DoubleObject<V, W> const& dob);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  template<class A, class B>
  DoubleObject<B> operator * (DoubleObject<A> lhs, DoubleObject<B> const& rhs) { return lhs *= rhs; }

  template<class T, class U>
  DoubleObject<T, U> operator * (double const& s, DoubleObject<T, U> rhs) { return rhs *= s; }

  template<class T, class U>
  DoubleObject<T, U> operator * (DoubleObject<T, U> lhs, double const& s) { return lhs *= s; }

  template<class T, class U>
  DoubleObject<T, U> operator / (double const& s, DoubleObject<T, U> rhs) { return rhs /= s; }

  template<class T, class U>
  DoubleObject<T, U> operator / (DoubleObject<T, U> lhs, double const& s) { return lhs /= s; }

  template<class T, class U>
  DoubleObject<T, U> operator * (DoubleObject<T, U> lhs, DoubleObject<T, U> const& rhs) { return lhs *= rhs; }

  template<class T, class U>
  DoubleObject<T, U> operator + (DoubleObject<T, U> lhs, DoubleObject<T, U> const& rhs) { return lhs += rhs; }

  template<class T, class U>
  DoubleObject<T, U> operator - (DoubleObject<T, U> lhs, DoubleObject<T, U> const& rhs) { return lhs -= rhs; }
  ///@}

  /**
   * @brief Method which prints the double object with cout <<.
   */
  template<class T, class U>
  std::ostream& operator << (std::ostream& os, DoubleObject<T, U> const& dob);
}
