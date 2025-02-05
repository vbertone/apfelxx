//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/set.h"
#include "apfel/operator.h"
#include "apfel/doubleobject.h"
#include "apfel/doubleoperator.h"
#include "apfel/distributionoperator.h"
#include "apfel/operatordistribution.h"
#include "apfel/messages.h"

#include <stdexcept>

namespace apfel
{
  //_________________________________________________________________________
  template<class T>
  Set<T>::Set(ConvolutionMap const& Map, std::map<int, T> const& in):
    _map(Map),
    _objects(in)
  {
  }

  //_________________________________________________________________________
  template<class T>
  Set<T>::Set(std::map<int, T> const& in):
    _map(ConvolutionMap{"DefaultConvolutionMap"}),
    _objects(in)
  {
    std::map<int, std::vector<ConvolutionMap::rule>> rules;
    for (auto const& o : in)
      rules.insert({o.first, {{o.first, o.first, 1}}});
    _map.SetRules(rules);
  }

  //_________________________________________________________________________
  template<class T>
  template<class V> Set<V> Set<T>::operator *= (Set<V> const& d) const
  {
    if (_map.GetName() != d.GetMap().GetName())
      throw std::runtime_error(error("Set::operator *=", "Convolution Map does not match"));

    std::map<int, V> mmap;
    for (auto item = _map.GetRules().begin(); item != _map.GetRules().end(); item++)
      {
        // If an element of the map with the same rules has already
        // been computed, retrieve it and use it.
        bool cycle = false;
        for (auto it = _map.GetRules().begin(); it != item; it++)
          if (it->second == item->second)
            {
              mmap.insert({item->first, mmap.at(it->first)});
              cycle = true;
              break;
            }
        if (cycle)
          continue;

        // Get set of distributions
        const auto& dist = d.GetObjects();

        // Start with the first object of the vector of rules. If it
        // does not exist continue.
        auto o = std::begin(item->second);
        if (dist.count(o->object) == 0)
          continue;
        V result = _objects.at(o->operand) * dist.at(o->object);

        // Multiply by the numerical coefficient only if it is
        // different from one.
        if(o->coefficient != 1)
          result *= o->coefficient;
        o++;

        // Continue with the following objects of the vector of rules
        for (auto end = std::end(item->second); o != end; o++)
          {
            // If the distribution does not exist skip it
            if (dist.count(o->object) == 0)
              continue;

            // Multiply by the numerical coefficient only if it is
            // different from one.
            if(o->coefficient == 0)
              continue;
            else if(o->coefficient != 1)
              result += o->coefficient * _objects.at(o->operand) * dist.at(o->object);
            else
              result += _objects.at(o->operand) * dist.at(o->object);
          }
        mmap.insert({item->first, result});
      }
    return Set<V> {d.GetMap(), mmap};
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
  Set<T>& Set<T>::operator *= (std::function<std::vector<double>(double const&)> f)
  {
    for (auto& v: _objects)
      v.second *= [=] (double const& x) -> double { return f(x)[v.first]; };

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
  Set<T>& Set<T>::operator *= (std::map<int, double> const& v)
  {
    for (auto& o: _objects)
      o.second *= v.at(o.first);

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
    // Initialize iterator on '_objects'
    auto it = _objects.begin();

    // Initialize 'CombObj' with the first object in '_objects'
    T CombObj = it->second;
    it++;

    // Continue with the following objects of the vector of rules
    for (auto end = _objects.end(); it != end; it++)
      CombObj += it->second;

    return CombObj;
  }

  //_________________________________________________________________________
  template<class T>
  T Set<T>::Combine(std::vector<double> const& weigths) const
  {
    // Check whether map of objects and vector of weights have the
    // same size.
    if (_objects.size() != weigths.size())
      throw std::runtime_error(error("Set::Combine", "Size of map of objects and vector of weights do not match"));

    // Initialize iterator on '_objects' and counter of the weights
    auto it = _objects.begin();
    int i;

    // In case the first weights are zero do not do the sum.
    for (i = 0; i < (int) weigths.size(); i++)
      if (weigths[i] != 0)
        break;
      else
        it++;

    // If the vector of weights is full of zeros just return the last
    // object multiplied by zero.
    if (i == (int) weigths.size())
      return 0 * (--it)->second;

    // Initialize 'CombObj' with the first object in '_objects'
    T CombObj = weigths[i++] * it->second;
    it++;

    // Continue with the following objects of the vector of rules
    for (auto end = _objects.end(); it != end; it++)
      {
        if (weigths[i] != 0)
          {
            if (weigths[i] == 1)
              CombObj += it->second;
            else
              CombObj += weigths[i] * it->second;
          }
        i++;
      }
    return CombObj;
  }

  //_________________________________________________________________________
  template<class U>
  std::ostream& operator << (std::ostream& os, Set<U> const& s)
  {
    os << "Set: " << &s << "\n";
    os << s.GetMap() << "\n";
    os << "Set of objects:\n";
    for (auto const& e : s.GetObjects())
      {
        os << "- Object index: " << e.first << "\n";
        os << "- Object:\n" << e.second << "\n";
      }
    return os;
  }

  // Specialisations
  template class Set<Distribution>;
  template class Set<Operator>;
  template class Set<DoubleObject<Distribution, Operator>>;
  template class Set<DoubleObject<Operator, Distribution>>;
  template class Set<DoubleDistribution>;
  template class Set<DoubleOperator>;
  template class Set<DistributionOperator>;
  template class Set<OperatorDistribution>;

  template Set<Distribution> Set<Operator>::operator *= (Set<Distribution> const&) const;
  template Set<Operator> Set<Operator>::operator *= (Set<Operator> const&) const;

  template std::ostream& operator << (std::ostream& os, Set<Distribution> const& s);
  template std::ostream& operator << (std::ostream& os, Set<Operator> const& s);
  template std::ostream& operator << (std::ostream& os, Set<DoubleObject<Distribution, Operator>> const& s);
  template std::ostream& operator << (std::ostream& os, Set<DoubleObject<Operator, Distribution>> const& s);
  template std::ostream& operator << (std::ostream& os, Set<DoubleDistribution> const& s);
  template std::ostream& operator << (std::ostream& os, Set<DoubleOperator> const& s);
  template std::ostream& operator << (std::ostream& os, Set<DistributionOperator> const& s);
  template std::ostream& operator << (std::ostream& os, Set<OperatorDistribution> const& s);
}
