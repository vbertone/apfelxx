//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

#include <functional>

using std::function;

namespace apfel
{
  /**
   * @brief The Set template class.
   *
   * This class allocates a Set of objects of type T
   * following the ConvolutionMap U. This class provides the
   * methods to deal with operations between different
   * types of objects T.
   *
   * Things to be improved:
   * - there is no real necessity of allocating U
   * - missing operators.
   */
  template<class T>
  class Set
  {        
  public:
    /**
     * @brief The Set class constructor.
     * @param the input map, in this case it makes a copy
     */
    Set(ConvolutionMap const& Map, map<int, T> const& in);

    /**
     * @brief operator *= product object
     * @param d left hand side object
     * @return a new object of type V and base U
     */
    template<class V> Set<V> operator *= (Set<V> const& d) const;

    // other operators
    Set<T>& operator *= (double const& s);                   //!< this *= scalar
    Set<T>& operator *= (function<double(double const&)> f); //!< this *= function of the integration variable (for distributions only)
    Set<T>& operator *= (vector<double> const& v);           //!< this *= vector of scalars
    Set<T>& operator /= (double const& s);                   //!< this /= scalar
    Set<T>& operator *= (Set<T> const& d);                   //!< this *= Set
    Set<T>& operator += (Set<T> const& d);                   //!< this += Set
    Set<T>& operator -= (Set<T> const& d);                   //!< this -= Set

    // Get methods
    T              const& at(int const& id) const { return _objects.at(id); }
    ConvolutionMap const& GetMap()          const { return _map; }
    map<int, T>    const& GetObjects()      const { return _objects; }

    // Method to sum all the objects of a given set.
    T Combine() const;

  private:
    ConvolutionMap _map;     //!< The shared pointer containing the convolution map
    map<int, T>    _objects; //!< The container for the map
  };

  /**
   * @brief operator *  and / definition
   * @param lhs the left object
   * @param rhs the right object
   * @return a Set of type B, C.
   */
  template<class A, class B>
  Set<B> operator * (Set<A> lhs, Set<B> const& rhs) { return lhs *= rhs; }

  // other operators
  template<class T>
  Set<T> operator * (double const& s, Set<T> rhs) { return rhs *= s; }

  template<class T>
  Set<T> operator * (Set<T> lhs, double const& s) { return lhs *= s; }

  template<class T>
  Set<T> operator * (function<double(double const&)> f, Set<T> rhs) { return rhs *= f; }

  template<class T>
  Set<T> operator * (Set<T> lhs, function<double(double const&)> f) { return lhs *= f; }

  template<class T>
  Set<T> operator * (vector<double> const& v, Set<T> rhs) { return rhs *= v; }

  template<class T>
  Set<T> operator * (Set<T> lhs, vector<double> const& v) { return lhs *= v; }

  template<class T>
  Set<T> operator / (int const& s, Set<T> rhs) { return rhs /= s; }

  template<class T>
  Set<T> operator / (Set<T> lhs, double const& s) { return lhs /= s; }

  template<class T>
  Set<T> operator * (Set<T> lhs, Set<T> const& rhs) { return lhs *= rhs; }

  template<class T>
  Set<T> operator + (Set<T> lhs, Set<T> const& rhs) { return lhs += rhs; }

  template<class T>
  Set<T> operator - (Set<T> lhs, Set<T> const& rhs) { return lhs -= rhs; }
}
