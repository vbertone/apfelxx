//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingbasisqcd.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/physicalbasisqcd.h"

namespace apfel
{
  //_________________________________________________________________________________
  MatchingBasisQCD::MatchingBasisQCD(int const& nf):
    ConvolutionMap{"MatchingBasisQCD_" + std::to_string(nf)}
  {
    // All valence-like distributions match multiplicatively through
    // M7.
    for (int k = 1; k <= 6; k++)
      _rules[2 * k] = {{M0, 2 * k, 1}, {M10, 2 * k, 1}};

    // Now we consider singlet like distributions
    const int nf1 = nf + 1;
    const int onf = 2 * nf1 - 1;

    // Gluon
    _rules[0] = {{M0, 0, 1}, {M1, 0, 1}, {M2, 1, 1./6.}, {M3, onf, - 1. / nf1}};
    for (int k = nf + 2; k <= 6; k++)
      _rules[0].push_back({M2, 2 * k - 1, 1. / k / ( k - 1 )});

    // Singlet
    _rules[1] = {{M0, 1, 1}, {M4, 0, 1}, {M5, 1, 1. / 6.}, {M6, onf, - 1. / nf1}};
    for (int k = nf + 2; k <= 6; k++)
      _rules[1].push_back({M5, 2 * k - 1, 1. / k / ( k - 1 )});

    // Light singlet-like distributions
    for (int k = 2; k <= nf; k++)
      _rules[2 * k - 1] = {{M0, 2 * k - 1, 1}, {M7, 2 * k - 1, 1}};

    // Heavy singlet-like distributions
    _rules[onf] = {{M0, onf, 1}, {M8, 0, - static_cast<double>(nf)}, {M9, 1, - nf / 6.}, {M7, 1, nf * nf1 / 6.}, {M6, onf, static_cast<double>(nf) / nf1}, {M7, onf, 1}, {M5, onf, 1. / nf1}, {M9, onf, - 1. / nf1}};
    for (int k = nf + 2; k <= 6; k++)
      {
        _rules[onf].push_back({M9, 2 * k - 1, - nf / (static_cast<double>(k * ( k - 1 )))});
        _rules[onf].push_back({M7, 2 * k - 1, nf * nf1 / (static_cast<double>(k * ( k - 1 )))});
      }

    // Super-heavy singlet-like distributions. They match exactly
    // like the singlet.
    for (int l = nf + 2; l <= 6; l++)
      {
        _rules[2 * l - 1] = {{M0, 2 * l - 1, 1}, {M4, 0, 1}, {M5, 1, 1./6.}, {M6, onf, - 1. / nf1}};
        for (int k = nf + 2; k <= 6; k++)
          _rules[2 * l - 1].push_back({M5, 2 * k - 1, 1. / k / ( k - 1 )});
      }

    // Flip sign of line corresponding to T3
    for (auto& t : _rules[3])
      t.coefficient *= -1;

    // Flip sign of columm corresponding to T3
    for (auto& r : _rules)
      for (auto& t : r.second)
        if (t.object == 3)
          t.coefficient *= -1;
  }

  //_________________________________________________________________________________
  MatchingOperatorBasisQCD::MatchingOperatorBasisQCD(int const& nf):
    ConvolutionMap{"MatchingOperatorBasisQCD_" + std::to_string(nf)}
  {
    // Allocate MatchingBasisQCD object to retrieve the splitting
    // matrix rules.
    const MatchingBasisQCD mb{nf};

    // Get matrix of coefficients
    const matrix<std::vector<double>> rc = mb.GetRuleMatrix();

    // Get matrix of operator indices
    const matrix<std::vector<int>> ri = mb.GetRuleIndices();

    // Now construct set of rules
    for (int i = 0; i < 13; i++)
      for (int j = 0; j < 13; j++)
        for (int k = 0; k < 13; k++)
          {
            if (rc(i, k).empty() || Gkj.count({k, j}) == 0)
              continue;

            for (int l = 0; l < (int) rc(i, k).size(); l++)
              _rules[Gkj.at({i, j})].push_back({ri(i, k)[l], Gkj.at({k, j}), rc(i, k)[l]});
          }
  }

  //_________________________________________________________________________________
  PhysicalMatchingBasisQCD::PhysicalMatchingBasisQCD(int const& nf):
    ConvolutionMap{"PhysicalMatchingBasisQCD_" + std::to_string(nf)}
  {
    // Leading-order matching
    for (int i = - 6 + 6; i <= 6 + 6; i++)
      _rules[i] = {{ONE, i, 1}};

    // Light minus-type distributions
    for (int i = - nf + 6; i <= - 1 + 6; i++)
      _rules[i].push_back({KLLm, i, 1});

    // Gluon
    _rules[GLUON].push_back({KGG, GLUON, 1});
    for (int j = 1 + 6; j <= nf + 6; j++)
      _rules[GLUON].push_back({KGL, j, 1});
    _rules[GLUON].push_back({KGH, nf + 1 + 6, 1});

    // Light plus-type distributions
    for (int i = 1 + 6; i <= nf + 6; i++)
      {
        _rules[i].push_back({KLG, 0 + 6, 1});
        _rules[i].push_back({KLL, i, 1});
        for (int j = 1 + 6; j <= nf + 6; j++)
          _rules[i].push_back({KLLP, j, 1});
      }

    // Heavy plus-type distribution
    _rules[nf + 1 + 6].push_back({KHG, 0 + 6, 1});
    _rules[nf + 1 + 6].push_back({KHH, nf + 1 + 6, 1});
    for (int j = 1 + 6; j <= nf + 6; j++)
      _rules[nf + 1 + 6].push_back({KHL, j, 1});
  }

  //_________________________________________________________________________________
  PhysicalMatchingOperatorBasisQCD::PhysicalMatchingOperatorBasisQCD(int const& nf):
    ConvolutionMap{"PhysicalMatchingOperatorBasisQCD_" + std::to_string(nf)}
  {
    // Allocate PhysicalMatchingBasisQCD object to retrieve the
    // splitting matrix rules.
    const PhysicalMatchingBasisQCD mb{nf};

    // Get matrix of coefficients
    const matrix<std::vector<double>> rc = mb.GetRuleMatrix();

    // Get matrix of operator indices
    const matrix<std::vector<int>> ri = mb.GetRuleIndices();

    // Now construct set of rules
    for (int i = 0; i < 13; i++)
      for (int j = 0; j < 13; j++)
        for (int k = 0; k < 13; k++)
          {
            if (rc(i, k).empty() || GkjPhys.count({k, j}) == 0)
              continue;

            for (int l = 0; l < (int) rc(i, k).size(); l++)
              _rules[GkjPhys.at({i, j})].push_back({ri(i, k)[l], GkjPhys.at({k, j}), rc(i, k)[l]});
          }
  }
}
