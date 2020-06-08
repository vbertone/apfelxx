//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/convolutionmap.h"

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
  DiagonalBasis::DiagonalBasis(int const& nf):
    ConvolutionMap{"DiagonalBasis_" + std::to_string(nf)}
  {
    for (int k = 0; k < nf; k++)
      _rules[k] = { {k, k, 1} };
  }
}
