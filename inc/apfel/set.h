//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/basismap.h>

#include <memory>
#include <unordered_map>
using std::unordered_map;
using std::shared_ptr;

namespace apfel
{
  /**
   * @brief The Set template class.
   *
   * This class allocates a Set of objects of type T
   * following the BasisMap U. This class provides the
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
    Set(BasisMap const& map, unordered_map<int, T> const& in);

    /**
     * @brief operator *= product object
     * @param d left hand side object
     * @return a new object of type V and base U
     */
    template<class V> Set<V> operator*=(Set<V> const& d) const;

    // other operators
    Set<T>& operator=(Set<T> const& d);  //!< this = Set<T>
    Set<T>& operator*=(double const& s); //!< this *= scalar
    Set<T>& operator*=(Set<T> const& d); //!< this *= Set
    Set<T>& operator+=(Set<T> const& d); //!< this += Set

    // Get methods
    T                     const& at(int const& id) const { return _objects.at(id); }
    BasisMap              const& GetMap()          const { return _map; }
    unordered_map<int, T> const& GetObjects()      const { return _objects; }

  private:
    BasisMap              const& _map;     //!< the shared pointer containin the flavor map
    unordered_map<int, T>        _objects; //!< The container for the unordered_map
  };

  /**
   * @brief operator * definition
   * @param lhs the left object
   * @param rhs the right object
   * @return a Set of type B, C.
   */
  template<class A, class B>
  inline Set<B> operator*(Set<A> lhs, Set<B> const& rhs) { return lhs *= rhs; }

  // other operators
  template<class T>
  inline Set<T> operator*(double const& s, Set<T> rhs) { return rhs *= s;}

  template<class T>
  inline Set<T> operator*(Set<T> lhs, double const& s) { return lhs *= s;}

  template<class T>
  inline Set<T> operator*(Set<T> lhs, Set<T> const& rhs) { return lhs *= rhs;}

  template<class T>
  inline Set<T> operator+(Set<T> lhs, Set<T> const& rhs) { return lhs += rhs;}

}
