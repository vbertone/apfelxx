//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"
#include "apfel/evolutionbasisqcd.h"

namespace apfel
{
  /**
   * @defgroup MatchBases Matching convolution maps
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for the matching at the heavy-quark
   * thresholds.
   * @ingroup ConvMap
   */
  ///@{
  /**
   * @brief The MatchingBasisQCD class is a derived of
   * ConvolutionMap specialised for the matching of distributions
   * using the QCD evolution basis and without assuming that intrinsic
   * heavy quark contributions vanish.
   */
  class MatchingBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {M1, M2, M3, M4, M5, M6, M7};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The MatchingBasisQCD constructor for the matching in the
     * QCD evolution basis with nf active flavours.
     @param nf: number of active flavours
     */
    MatchingBasisQCD(int const& nf):
      ConvolutionMap{"MatchingBasisQCD_" + std::to_string(nf)}
    {
      // All valence-like distributions match multiplicatively through
      // M7.
      for (int k = 1; k <= 6; k++)
        _rules[2 * k] = { {M7, 2 * k, 1} };

      // Now we consider singlet like distributions.
      const int nf1 = nf + 1;
      const int onf = 2 * nf1 - 1;

      // Gluon
      _rules[0] = { {M1, 0, 1}, {M2, 1, 1./6.}, {M3, onf, - 1. / nf1} };
      for (int k = nf + 2; k <= 6; k++)
        _rules[0].push_back({M2, 2 * k - 1, 1. / k / ( k - 1 )});

      // Singlet
      _rules[1] = { {M4, 0, 1}, {M5, 1, 1./6.}, {M6, onf, - 1. / nf1} };
      for (int k = nf + 2; k <= 6; k++)
        _rules[1].push_back({M5, 2 * k - 1, 1. / k / ( k - 1 )});

      // Light singlet-like distributions
      for (int k = 2; k <= nf; k++)
        _rules[2 * k - 1] = { {M7, 2 * k - 1, 1} };

      // Heavy singlet-like distributions
      _rules[onf] = { {M4, 0, - static_cast<double>(nf)}, {M5, 1, - nf / 6.}, {M7, 1, nf * nf1 / 6.}, {M6, onf, static_cast<double>(nf) / nf1}, {M7, onf, 1} };
      for (int k = nf + 2; k <= 6; k++)
        {
          _rules[onf].push_back({M5, 2 * k - 1, - nf / (static_cast<double>(k * ( k - 1 )))});
          _rules[onf].push_back({M7, 2 * k - 1, nf * nf1 / (static_cast<double>(k * ( k - 1 )))});
        }

      // Super-heavy singlet-like distributions. They match exactly
      // like the singlet.
      for (int l = nf + 2; l <= 6; l++)
        {
          _rules[2 * l - 1] = { {M4, 0, 1}, {M5, 1, 1./6.}, {M6, onf, - 1. / nf1} };
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
    };
  };

  /**
   * @brief The MatchingOperatorBasisQCD class is a derived of
   * ConvolutionMap specialised for the matching of the evolution of
   * operators at the heavy-quark thresholds using the QCD evolution
   * basis.
   */
  class MatchingOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The MatchingOperatorBasisQCD constructor
     * @param nf: number of active flavours
     */
    MatchingOperatorBasisQCD(int const& nf):
      ConvolutionMap{"MatchingOperatorBasisQCD_" + std::to_string(nf)}
    {
      // Allocate MatchingBasisQCD object to retrieve the splitting
      // matrix rules
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
    };
  };
  ///@}
}
