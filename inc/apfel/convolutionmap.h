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
     * @brief ConvolutionMap constructor
     * @param name: name of the map
     */
    ConvolutionMap(std::string const& name);

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
    std::map<int,std::vector<rule>> const& GetRules() const { return _rules; }

    /**
     * @brief Retrieve the full set of rules for the multiplications
     * in the form of a matrix.
     * @return The multiplication rule matrix
     */
    matrix<double> const GetRuleMatrix() const
    {
      matrix<double> m{_rules.size(), _rules.size()};
      for (auto const& r : _rules)
        for (auto const& e : r.second)
          m(r.first, e.object) = e.coefficient;

      return m;
    }
    ///@}
  protected:
    std::map<int,std::vector<rule>> _rules; //!< the map container
    std::string                     _name;  //!< the name of the derived class
  };
}
