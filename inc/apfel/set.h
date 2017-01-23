//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <memory>
#include <unordered_map>
#include <apfel/basismap.h>
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
  template<class T, class U>
  class Set
  {        
  public:
    /**
     * @brief The Set class constructor.
     * @param the input map.
     */
    Set(unordered_map<int, T> const& in):
      _objects(in)
    {
      _map = unique_ptr<U>(new U{});
    }

    /**
     * @brief operator *= product object
     * @param d left hand side object
     * @return a new object of type V and base U
     */
    template<class V>
    Set<V,U> operator*=(Set<V,U> const& d) const
    {
      unordered_map<int,V> mmap;
      for (auto const& item: _map->GetRules())
        {
          auto o = std::begin(item.second);
          V result = (*o).operation*_objects.at((*o).splitting)*d.GetObjects().at((*o).distribution);
          o++;
          for (auto end = std::end(item.second); o != end; o++)
            result += (*o).operation*_objects.at((*o).splitting)*d.GetObjects().at((*o).distribution);
          mmap.insert({item.first,result});
        }
      return Set<V,U>{mmap};
    }

    // Get methods
    unordered_map<int, T> const& GetObjects() const { return _objects; }

  private:
    unordered_map<int, T> const& _objects; //!< The container for the unordered_map
    shared_ptr<U> _map;                    //!< the shared pointer containin the flavor map
  };

  /**
   * @brief operator * definition
   * @param lhs the left object
   * @param rhs the right object
   * @return a Set of type B, C.
   */
  template<class A, class B, class C>
  Set<B,C> operator*(Set<A,C> lhs, Set<B,C> const& rhs)
  {
    return lhs *= rhs;
  }


}
