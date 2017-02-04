//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <unordered_map>
#include <string>
#include <vector>
using std::unordered_map;
using std::string;
using std::vector;

namespace apfel
{
  /**
   * @brief The BasisMap class provides a common set of tools
   * to autodiscover the combination rules between sets of
   * operators and distributions.
   */
  class BasisMap
  {
  public:
    /**
     * @brief BasisMap
     * @param name
     */
    BasisMap(string const& name);

    /**
     * @brief The rule struct
     */
    struct rule
    {
      int    operand;
      int    object;
      double coefficient;
    };

    // Get methods
    string                          const& GetName()  const { return _name; }
    unordered_map<int,vector<rule>> const& GetRules() const { return _rules; }

  protected:
    unordered_map<int,vector<rule>> _rules; //!< the map container
    string _name;                           //!< the name of the derived class
  };

}
