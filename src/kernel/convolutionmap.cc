//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/convolutionmap.h"

#include <sstream>
#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  ConvolutionMap::ConvolutionMap(std::string const& name):
    _name(name)
  {
  }

  //_________________________________________________________________________________
  matrix<std::vector<double>> const ConvolutionMap::GetRuleMatrix() const
  {
    matrix<std::vector<double>> m;
    m.resize(_rules.size(), _rules.size(), {});
    for (auto const& r : _rules)
      for (auto const& e : r.second)
        m(r.first, e.object).push_back(e.coefficient);

    return m;
  }

  //_________________________________________________________________________________
  matrix<std::vector<int>> const ConvolutionMap::GetRuleIndices() const
  {
    matrix<std::vector<int>> m{};
    m.resize(_rules.size(), _rules.size(), {});
    for (auto const& r : _rules)
      for (auto const& e : r.second)
        m(r.first, e.object).push_back(e.operand);

    return m;
  }

  //_________________________________________________________________________________
  std::ostream& operator << (std::ostream& os, ConvolutionMap const& cm)
  {
    os << "ConvolutionMap: " << &cm << "\n";
    os << "Name: " << cm.GetName() << "\n";
    os << "Operand-index matrix:\n";
    const matrix<std::vector<int>> ri = cm.GetRuleIndices();
    for (int i = 0; i < (int) ri.size(0); i++)
      {
        os << i << ":\t";
        for (int j = 0; j < (int) ri.size(1); j++)
          {
            os << "{";
            for (auto const& e : ri(i, j))
              os << e << ", ";
            if (!ri(i, j).empty())
              os << "\b\b";
            os << "} ";
          }
        os << "\n";
      }
    os << "Coefficient matrix:\n";
    const matrix<std::vector<double>> rc = cm.GetRuleMatrix();
    const std::ostringstream default_format;
    os << std::scientific;
    os.precision(1);
    for (int i = 0; i < (int) rc.size(0); i++)
      {
        os << i << ":\t";
        for (int j = 0; j < (int) rc.size(1); j++)
          {
            os << "{";
            for (auto const& e : rc(i, j))
              os << e << ", ";
            if (!rc(i, j).empty())
              os << "\b\b";
            os << "} ";
          }
        if (i != (int) rc.size(0) - 1)
          os << "\n";
      }
    os.copyfmt(default_format);
    return os;
  }

  //_________________________________________________________________________________
  ConvolutionMap::rule operator * (double const& s, ConvolutionMap::rule rhs)
  {
    return rhs *= s;
  }

  //_________________________________________________________________________________
  ConvolutionMap::rule operator * (ConvolutionMap::rule lhs, double const& s)
  {
    return lhs *= s;
  }

  //_________________________________________________________________________________
  ConvolutionMap::rule operator / (ConvolutionMap::rule lhs, double const& s)
  {
    return lhs /= s;
  }

  //_________________________________________________________________________________
  std::vector<ConvolutionMap::rule> operator + (std::vector<ConvolutionMap::rule> lhs, std::vector<ConvolutionMap::rule> rhs)
  {
    // Start with the l.h.s.
    std::vector<ConvolutionMap::rule> sum = lhs;

    // Now run over the entries of the r.h.s.
    for (ConvolutionMap::rule const& r : rhs)
      {
        // Find iterator of r in sum
        const auto it = std::find_if(sum.begin(), sum.end(), [&r] (ConvolutionMap::rule const& r1) -> bool { return (r.operand == r1.operand && r.object == r1.object ? true : false); });

        // If there is no rule in sum with the same operand and object
        // push the new rule back into the sum ...
        if (it == sum.end())
          sum.push_back(r);
        // ... otherwise sum the coefficient to the existing one
        else
          (*it).coefficient += r.coefficient;
      }
    return sum;
  }

  //_________________________________________________________________________________
  std::vector<ConvolutionMap::rule> operator - (std::vector<ConvolutionMap::rule> lhs, std::vector<ConvolutionMap::rule> rhs)
  {
    // Start with the l.h.s.
    std::vector<ConvolutionMap::rule> diff = lhs;

    // Now run over the entries of the r.h.s.
    for (ConvolutionMap::rule const& r : rhs)
      {
        // Find iterator of r in diff
        const auto it = std::find_if(diff.begin(), diff.end(), [&r] (ConvolutionMap::rule const& r1) -> bool { return (r.operand == r1.operand && r.object == r1.object ? true : false); });

        // If there is no rule in diff with the same operand and object
        // push the new rule back into the diff with a minus sign ...
        if (it == diff.end())
          diff.push_back((-1)*r);
        // ... otherwise subtract the coefficient from the existing one
        else
          (*it).coefficient -= r.coefficient;
      }
    return diff;
  }

  //_________________________________________________________________________________
  std::vector<ConvolutionMap::rule> operator * (double const& s, std::vector<ConvolutionMap::rule> rhs)
  {
    for (ConvolutionMap::rule& r : rhs)
      r *= s;
    return rhs;
  }

  //_________________________________________________________________________________
  std::vector<ConvolutionMap::rule> operator * (std::vector<ConvolutionMap::rule> lhs, double const& s)
  {
    for (ConvolutionMap::rule& r : lhs)
      r *= s;
    return lhs;
  }

  //_________________________________________________________________________________
  std::vector<ConvolutionMap::rule> operator / (std::vector<ConvolutionMap::rule> lhs, double const& s)
  {
    for (ConvolutionMap::rule& r : lhs)
      r /= s;
    return lhs;
  }
}
