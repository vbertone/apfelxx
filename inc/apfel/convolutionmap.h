//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

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
    std::string const& GetName()  const { return _name; }

    /**
     * @brief Retrieve the full set of rules for the multiplications.
     * @return The multiplication rules
     */
    std::map<int,std::vector<rule>> const& GetRules() const { return _rules; }
    ///@}
  protected:
    std::map<int,std::vector<rule>> _rules; //!< the map container
    std::string                     _name;  //!< the name of the derived class
  };
}
