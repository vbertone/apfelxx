//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/convolutionmap.h"

#include <iostream>
#include <sstream>

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
    os << "Operator-index matrix:\n";
    const matrix<std::vector<int>> ri = cm.GetRuleIndices();
    for (int i = 0; i < (int) ri.size(0); i++)
      {
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
  DiagonalBasis::DiagonalBasis(int const& nf):
    ConvolutionMap{"DiagonalBasis_" + std::to_string(nf)}
  {
    for (int k = 0; k < nf; k++)
      _rules[k] = { {k, k, 1} };
  }
}
