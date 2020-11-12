//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

#include <functional>

namespace apfel
{
  /**
   * @brief The Set template class allocates a collection of objects
   * of type T along the ConvolutionMap and provides the methods to
   * perform operations between different types of objects T.
   */
  template<class T>
  class Set
  {
  public:
    /**
     * @brief The Set constructor.
     * @param Map: the convolution map (defualt: empty ConvolutionMap)
     * @param in: a map of objects of type T (defualt: empty map of objects)
     */
    Set(ConvolutionMap const& Map = ConvolutionMap{"UNDEFINED"}, std::map<int, T> const& in = std::map<int, T>{});

    /**
     * @name Binary operators
     */
    ///@{
    /**
     * @brief operator *= product object
     * @param d: the left hand side object of type V
     * @return a new object of type V
     */
    template<class V> Set<V> operator *= (Set<V> const& d) const;

    Set<T>& operator *= (double const& s);                             //!< this *= scalar
    Set<T>& operator *= (std::function<double(double const&)> f);      //!< this *= function of the integration variable (for distributions only)
    Set<T>& operator *= (std::vector<double> const& v);                //!< this *= vector of scalars
    Set<T>& operator /= (double const& s);                             //!< this /= scalar
    Set<T>& operator += (Set<T> const& d);                             //!< this += Set
    Set<T>& operator -= (Set<T> const& d);                             //!< this -= Set
    ///@}

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief This returns object with ID "id" in the map.
     */
    T const& at(int const& id) const { return _objects.at(id); }
    /**
     * @brief This returns the convolution map.
     */
    ConvolutionMap const& GetMap() const { return _map; }
    /**
     * @brief This returns the full map of objects.
     */
    std::map<int, T> const& GetObjects() const { return _objects; }
    ///@}

    /**
     * @brief This function (re)sets the convolution map.
     */
    void SetMap(ConvolutionMap const& map) { _map = map; }

    /**
     * @brief This function (re)sets the map of objects.
     */
    void SetObjects(std::map<int, T> const& objects) { _objects = objects; }

    /**
     * @brief This function sums up all the objects of the set into
     * one.
     */
    T Combine() const;

    /**
     * @brief This function sums up all the objects of the set into
     * one using the components of input vactor for as weights.
     */
    T Combine(std::vector<double> const& weigths) const;

  private:
    ConvolutionMap   _map;     //!< The shared pointer containing the convolution map
    std::map<int, T> _objects; //!< The container for the map

    template<class U>
    friend std::ostream& operator << (std::ostream& os, Set<U> const& s);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  template<class A, class B>
  Set<B> operator * (Set<A> lhs, Set<B> const& rhs) { return lhs *= rhs; }

  // other operators
  template<class T>
  Set<T> operator * (double const& s, Set<T> rhs) { return rhs *= s; }

  template<class T>
  Set<T> operator * (Set<T> lhs, double const& s) { return lhs *= s; }

  template<class T>
  Set<T> operator * (std::function<double(double const&)> f, Set<T> rhs) { return rhs *= f; }

  template<class T>
  Set<T> operator * (Set<T> lhs, std::function<double(double const&)> f) { return lhs *= f; }

  template<class T>
  Set<T> operator * (std::vector<double> const& v, Set<T> rhs) { return rhs *= v; }

  template<class T>
  Set<T> operator * (Set<T> lhs, std::vector<double> const& v) { return lhs *= v; }

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
  ///@}
}
