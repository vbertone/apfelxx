//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matrix.h"

#include <map>
#include <string>
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
      rule& operator *= (double const& s)
      {
        coefficient *= s;
        return *this;
      }
      rule& operator /= (double const& s)
      {
        coefficient /= s;
        return *this;
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
   * @name Ternary operators on ConvolutionMap::rule objects
   */
  ///@{
  ConvolutionMap::rule operator * (double const& s, ConvolutionMap::rule rhs);
  ConvolutionMap::rule operator * (ConvolutionMap::rule lhs, double const& s);
  ConvolutionMap::rule operator / (ConvolutionMap::rule lhs, double const& s);
  std::vector<ConvolutionMap::rule> operator + (std::vector<ConvolutionMap::rule> lhs, std::vector<ConvolutionMap::rule> rhs);
  std::vector<ConvolutionMap::rule> operator - (std::vector<ConvolutionMap::rule> lhs, std::vector<ConvolutionMap::rule> rhs);
  std::vector<ConvolutionMap::rule> operator * (double const& s, std::vector<ConvolutionMap::rule> lhs);
  std::vector<ConvolutionMap::rule> operator * (std::vector<ConvolutionMap::rule> lhs, double const& s);
  std::vector<ConvolutionMap::rule> operator / (std::vector<ConvolutionMap::rule> lhs, double const& s);
  ///@}
}
