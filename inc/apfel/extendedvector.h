//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>

namespace apfel
{
  /**
   * @brief Class that extends vectors to negative indices
   */
  template<class T>
  class ExtendedVector
  {
  public:
    /**
     * @brief The ExtendedVector constructor.
     * @param size: the size of the container (default: 0)
     * @param value: initialisation value (default: 0)
     * @param imin: the lowest allowed index (default: 0)
     */
    ExtendedVector(int const& size = 0, T const& value = 0, int const& imin = 0): _imin(imin), _vector(size, value) {}

    /**
     * @brief Returns the value at the possibly negative given index
     * @param index: the position index
     * @return the value at the possibly negative given index
     * @note Settable version
     */
    T& operator [] (int const& index) { return _vector[index - _imin]; }

    /**
     * @brief Returns the value at the possibly negative given index
     * @param index: the position index
     * @return the value at the possibly negative given index
     * @note Non-settable version
     */
    const T& operator [] (int const& index) const { return _vector[index - _imin]; }

    /**
     * @brief Returns the lower bound
     * @return the lower bound
     */
    int min() const { return _imin; };

    /**
     * @brief Returns the upper bound
     * @return the upper bound
     */
    int max() const { return _vector.size() + _imin; };

    /**
     * @brief Returns the size of the vector
     */
    size_t size() const { return _vector.size(); }

    /**
     * @brief Resizes the continer
     * @param size: the new size
     * @param value: the value used to fill in the additional (if any) slots
     * @param imin: the lowest allowed index (default: 0)
     */
    void resize(int const& size, T const& value = 0, int const& imin = 0) { _imin = imin; _vector.resize(size, value); }

    /**
     * @brief Non-constant begin iterator
     */
    typename std::vector<T>::iterator begin() { return _vector.begin(); }

    /**
     * @brief Constant begin iterator
     */
    typename std::vector<T>::const_iterator begin() const { return _vector.begin(); }

    /**
     * @brief Non-constant end iterator
     */
    typename std::vector<T>::iterator end() { return _vector.end(); }

    /**
     * @brief Constant end iterator
     */
    typename std::vector<T>::const_iterator end() const { return _vector.end(); }

  private:
    int            _imin;    //!< The lower bound
    std::vector<T> _vector;  //!< The container
  };
}
