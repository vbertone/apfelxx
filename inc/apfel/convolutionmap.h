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
   * @brief The ConvolutionMap class provides a common set of tools to
   * autodiscover the combination rules between sets of operators and
   * distributions.
   */
  class ConvolutionMap
  {
  public:
    /**
     * @brief ConvolutionMap
     * @param name
     */
    ConvolutionMap(string const& name);

    /**
     * @brief The rule struct
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

    // Get methods
    string                const& GetName()  const { return _name; }
    map<int,vector<rule>> const& GetRules() const { return _rules; }

  protected:
    map<int,vector<rule>> _rules; //!< the map container
    string                _name;  //!< the name of the derived class
  };
}
