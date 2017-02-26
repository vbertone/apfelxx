//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/tools.h"

#include <algorithm>
#include <list>
#include <mutex>
#include <thread>
#include <unordered_map>

using std::list;
using std::mutex;
using std::unordered_map;
using std::lock_guard;
using std::pair;

namespace std
{
  namespace {

    //! Code from boost
    //! Reciprocal of the golden ratio helps spread entropy and handles duplicates.
    //! See Mike Seymour in magic-numbers-in-boosthash-combine:
    //! http://stackoverflow.com/questions/4948780
    template <class T> inline void hash_combine(std::size_t& seed, T const& v)
    {
      seed ^= hash<T>()(v) + 0x9e3779b9 + (seed<<6) + (seed>>2);
    }

    //! Recursive template code derived from Matthieu M.
    template <class Tuple, size_t Index = std::tuple_size<Tuple>::value - 1>
    struct HashValueImpl
    {        
      static void apply(size_t& seed, Tuple const& tuple)
      {
        HashValueImpl<Tuple, Index-1>::apply(seed, tuple);
        hash_combine(seed, get<Index>(tuple));
      }
    };

    template <class Tuple>
    struct HashValueImpl<Tuple,0>
    {
      static void apply(size_t& seed, Tuple const& tuple)
      {
        hash_combine(seed, get<0>(tuple));
      }
    };
  }

  //! Hashing of tuples.
  //! From http://stackoverflow.com/questions/7110301/generic-hash-for-tuples-in-unordered-map-unordered-set
  template <typename ... TT>
  struct hash<std::tuple<TT...>>
  {
    size_t operator()(std::tuple<TT...> const& tt) const
    {
      size_t seed = 0;
      HashValueImpl<std::tuple<TT...> >::apply(seed, tt);
      return seed;
    }
  };
}

namespace apfel {

  /**
   * @brief The Cache class
   *
   * This class provide a simple LRU implementation based on unordered_map.
   */
  template <class Key, class Value>
  class Cache {
   public:

    /**
     * @brief Cache
     * @param maxSize
     * @param elasticity
     */
    Cache(size_t maxSize = 100, size_t elasticity = 10):
      _maxsize(maxSize), _maxallowed(maxSize+elasticity) {}

    //! delete copy constructor
    Cache(Cache const&) = delete;

    /**
     * @brief size
     * @return
     */
    size_t size() const
    {
      lock_guard<mutex> g(lock_);
      return _cache.size();
    }

    /**
     * @brief empty
     * @return
     */
    bool empty() const
    {
      lock_guard<mutex> g(lock_);
      return _cache.empty();
    }

    /**
     * @brief clear cache (map and keys)
     */
    void clear()
    {
      lock_guard<mutex> g(lock_);
      _cache.clear();
      _keys.clear();
    }

    /**
     * @brief insert value for key
     * @param k the key
     * @param v the value
     */
    void insert(const Key& k, const Value& v)
    {
      lock_guard<mutex> g(lock_);
      const auto iter = _cache.find(k);
      if (iter != _cache.end())
        {
          iter->second->second = v;
          _keys.splice(_keys.begin(), _keys, iter->second);
          return;
        }
      _keys.emplace_front(k, v);
      _cache[k] = _keys.begin();
      prune();
    }

    /**
     * @brief Try and get key if available
     * @param kIn the key
     * @param vOut the output value if available
     * @return
     */
    bool tryGet(const Key& kIn, Value& vOut)
    {
      lock_guard<mutex> g(lock_);
      const auto iter = _cache.find(kIn);
      if (iter == _cache.end()) return false;
      _keys.splice(_keys.begin(), _keys, iter->second);
      vOut = iter->second->second;
      return true;
    }

    /**
     * @brief get value for key
     * @param k the key
     * @return the value
     */
    const Value& get(const Key& k)
    {
      lock_guard<mutex> g(lock_);
      const auto iter = _cache.find(k);
      if (iter == _cache.end()) {
        throw runtime_exception("Cache::Value","key not found.");
      }
      _keys.splice(_keys.begin(), _keys, iter->second);
      return iter->second->second;
    }

    /**
     * @brief remove key
     * @param k the key
     */
    bool remove(const Key& k)
    {
      lock_guard<mutex> g(lock_);
      auto iter = _cache.find(k);
      if (iter == _cache.end()) return false;
      _keys.erase(iter->second);
      _cache.erase(iter);
      return true;
    }

    /**
     * @brief check if cache contains key
     * @param k the key
     * @return true if key is there, false otherwise
     */
    bool contains(const Key& k)
    {
      lock_guard<mutex> g(lock_);
      return _cache.find(k) != _cache.end();
    }

    // Getters
    size_t getMaxSize()        const { return _maxsize; }    //!< return the hard limit
    size_t getMaxAllowedSize() const { return _maxallowed; } //!< return the soft limit

  private:
    /**
     * @brief prune objects when size of cache saturates
     * @return
     */
    void prune()
    {
      if (_maxsize == 0 || _cache.size() < _maxallowed) return;

      size_t count = 0;
      while (_cache.size() > _maxsize)
        {
          _cache.erase(_keys.back().first);
          _keys.pop_back();
          count++;
        }
    }

    mutable mutex lock_;
    unordered_map<Key, typename list<pair<Key, Value>>::iterator> _cache;
    list<pair<Key, Value>> _keys;  //!< store the key, value map
    size_t _maxsize, _maxallowed;
  };

}

