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
    ///@}
  protected:
    std::map<int, std::vector<rule>> _rules; //!< the map container
    std::string                      _name;  //!< the name of the derived class
  };

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
     */
    DiagonalBasis(int const& nf);
  };
}
