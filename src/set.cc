//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/set.h"
#include "apfel/distribution.h"
#include "apfel/operator.h"
#include "apfel/tools.h"

namespace apfel {

  //_________________________________________________________________________
  template<class T>
  Set<T>::Set(BasisMap const& map, const unordered_map<int,T> &in):
    _map(map),
    _objects(in)
  {
  }

  //_________________________________________________________________________
  template<class T>
  template<class V> Set<V> Set<T>::operator*=(Set<V> const& d) const
  {
    if (&_map != &d.GetMap())
      throw runtime_exception("Set::operator*=", "Basis Map grid does not match");

    unordered_map<int,V> mmap;
    for (auto const& item: _map.GetRules())
      {
        auto o = std::begin(item.second);
        V result = (*o).operation*_objects.at((*o).splitting)*d.GetObjects().at((*o).distribution);
        o++;
        for (auto end = std::end(item.second); o != end; o++)
          result += (*o).operation*_objects.at((*o).splitting)*d.GetObjects().at((*o).distribution);
        mmap.insert({item.first,result});
      }

    return Set<V>{d.GetMap(),mmap};
  }

  template class Set<Distribution>;
  template class Set<Operator>;
  template Set<Distribution> Set<Operator>::operator*=(Set<Distribution> const&) const;
}
