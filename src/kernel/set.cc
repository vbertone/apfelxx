//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/set.h"
#include "apfel/operator.h"
#include "apfel/messages.h"

#include <stdexcept>

namespace apfel {
  //_________________________________________________________________________
  template<class T>
  Set<T>::Set(ConvolutionMap const& Map, std::map<int,T> const& in):
    _map(Map),
    _objects(in)
  {
  }

  //_________________________________________________________________________
  template<class T>
  template<class V> Set<V> Set<T>::operator *= (Set<V> const& d) const
  {
    if (_map.GetName() != d.GetMap().GetName())
      throw std::runtime_error(error("Set::operator *=", "Convolution Map does not match (1)"));

    std::map<int,V> mmap;
    for (auto item = _map.GetRules().begin(); item != _map.GetRules().end(); item++)
      {
	// If an element of the map with the same rules has already
	// been computed, retrieve it and use it.
	bool cycle = false;
	for (auto it = _map.GetRules().begin(); it != item; it++)
	  if (it->second == item->second)
	    {
	      mmap.insert({item->first,mmap.at(it->first)});
	      cycle = true;
	      break;
	    }
	if (cycle)
	  continue;

	// Get set of distributions.
	const auto& dist = d.GetObjects();

	// Start with the first object of the vector or rules.
	// If it does not exist, continue.
        auto o = std::begin(item->second);
	if (dist.count((*o).object) == 0)
	  continue;
        V result = _objects.at((*o).operand) * dist.at((*o).object);

	// Multiply by the numerical coefficient only if it is
	// different from one.
	if((*o).coefficient != 1)
	  result *= (*o).coefficient;
        o++;

	// Continue with the following objects of the vector of rules.
        for (auto end = std::end(item->second); o != end; o++)
	  {
	    // If the distribution does not exist skip it.
	    if (dist.count((*o).object) == 0)
	      continue;

	    // Multiply by the numerical coefficient only if it is
	    // different from one.
	    if((*o).coefficient == 0)
	      continue;
	    else if((*o).coefficient != 1)
	      result += (*o).coefficient * _objects.at((*o).operand) * dist.at((*o).object);
	    else
	      result += _objects.at((*o).operand) * dist.at((*o).object);
	  }
        mmap.insert({item->first,result});
      }

    return Set<V>{d.GetMap(),mmap};
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator *= (double const& s)
  {
    for (auto& v: _objects)
      v.second *= s;
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator *= (std::function<double(double const&)> f)
  {
    for (auto& v: _objects)
      v.second *= f;
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator *= (std::vector<double> const& v)
  {
    for (auto& o: _objects)
      o.second *= v[o.first];
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator /= (double const& s)
  {
    const double r = 1. / s;
    for (auto& v: _objects)
      v.second *= r;
    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator += (Set<T> const& d)
  {
    if (_map.GetName() != d.GetMap().GetName())
      throw std::runtime_error(error("Set::operator +=", "Convolution Map does not match"));

    for (auto& v: _objects)
      v.second += d.at(v.first);

    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>& Set<T>::operator -= (Set<T> const& d)
  {
    if (_map.GetName() != d.GetMap().GetName())
      throw std::runtime_error(error("Set::operator -=", "Convolution Map does not match"));

    for (auto& v: _objects)
      v.second -= d.at(v.first);

    return *this;
  }

  //_________________________________________________________________________
  template<class T>
  T Set<T>::Combine() const
  {
    // Initialize iterator on '_objects'.
    auto it = _objects.begin();

    // Initialize 'CombObj' with the first object in '_objects'.
    T CombObj = it->second;
    it++;

    // Continue with the following objects of the vector of rules.
    for (auto end = _objects.end(); it != end; it++)
      CombObj += it->second;

    return CombObj;
  }

  // Specialisations.
  template class Set<Distribution>;
  template class Set<Operator>;
  template Set<Distribution> Set<Operator>::operator *= (Set<Distribution> const&) const;
  template Set<Operator> Set<Operator>::operator *= (Set<Operator> const&) const;
}
