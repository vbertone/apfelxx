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

using std::map;
using std::string;
using std::vector;

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
     * @brief ConvolutionMap default constructor
     * @param name: name of the map
     */
    ConvolutionMap(string const& name);

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
     * @brief Retrieve the name of the map.
     * @return The name of the map
     */
    string                const& GetName()  const { return _name; }

    /**
     * @brief Retrieve the full set of rules for the multiplications.
     * @return The multiplication rules
     */
    map<int,vector<rule>> const& GetRules() const { return _rules; }

  protected:
    map<int,vector<rule>> _rules; //!< the map container
    string                _name;  //!< the name of the derived class
  };
}
