//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matrix.h"

#include <map>
#include <string>
#include <vector>
#include <iostream>

namespace apfel
{
  /**
   * @brief The ConvolutionMap class encapsulates the set of rules to
   * multiply a sets of operators with a set of distributions.
   */
  class ConvolutionMap
  {
  public:
    /**
     * @brief This structure contains the attribute of a single rule.
     */
    struct rule
    {
      int    operand;
      int    object;
      double coefficient;
      bool   operator == (rule const& r) const
      {
        if (r.operand     != operand)     return false;
        if (r.object      != object)      return false;
        if (r.coefficient != coefficient) return false;
        return true;
      }
    };

    /**
     * @brief ConvolutionMap constructor
     * @param name: name of the map
     */
    ConvolutionMap(std::string const& name);

    /**
     * @name Setters
     */
    ///@{
    /**
     * @brief Set the rule of the convolution map.
     * @param rules: the input set of rules
     */
    void SetRules(std::map<int, std::vector<rule>> const& rules)  { _rules = rules; }
    ///@}

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Retrieve the name of the map.
     * @return The name of the map
     */
    std::string const& GetName() const { return _name; }

    /**
     * @brief Retrieve the full set of rules for the multiplications.
     * @return The multiplication rules
     */
    std::map<int, std::vector<rule>> const& GetRules() const { return _rules; }

    /**
     * @brief Retrieve the full set of rules for the multiplications
     * in the form of a matrix.
     */
    matrix<std::vector<double>> const GetRuleMatrix() const;

    /**
     * @brief Retrieve the operand indices of the full set of rules
     * for the multiplications in the form of a matrix. Elements set
     * to -1 correspond to empty slots.
     */
    matrix<std::vector<int>> const GetRuleIndices() const;

    /**
     * @brief Print the Operator object
     */
    void Print() const { std::cout << *this << std::endl; }
    ///@}
  protected:
    std::map<int, std::vector<rule>> _rules; //!< the map container
    std::string                      _name;  //!< the name of the derived class

    friend std::ostream& operator << (std::ostream& os, ConvolutionMap const& cm);
  };

  /**
   * @brief Method which prints ConvolutionMap with cout <<.
   */
  std::ostream& operator << (std::ostream& os, ConvolutionMap const& cm);

  /**
   * @brief The DiagonalBasis class is the simplest derivation of
   * ConvolutionMap meant to essentially perform a scalar product of
   * two sets of objects.
   */
  class DiagonalBasis: public ConvolutionMap
  {
  public:
    /**
     * @brief The DiagonalBasis constructor
     * @param nf: number of elements
     * @param offset: starting index for the enumeration on the distributions (default: 0)
     */
    DiagonalBasis(int const& nf, int const& offset = 0);
  };
}
