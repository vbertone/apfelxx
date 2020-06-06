//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

namespace apfel
{
  /**
   * @defgroup MatchBases Matching convolution maps Collection of
   * derived classes from ConvolutionMap that implement the
   * convolution map for the matching at the heavy-quark thresholds.
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
    enum Operand: int {M0, M1, M2, M3, M4, M5, M6, M7};
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
      // M0 and M7.
      for (int k = 1; k <= 6; k++)
        _rules[2 * k] = { {M0, 2 * k, 1}, {M7, 2 * k, 1} };

      // Now we consider singlet like distributions.
      const int nf1 = nf + 1;
      const int onf = 2 * nf1 - 1;

      // Gluon
      _rules[0] = { {M0, 0, 1}, {M1, 0, 1}, {M2, 1, 1./6.}, {M3, onf, - 1. / nf1} };
      for (int k = nf + 2; k <= 6; k++)
        _rules[0].push_back({M2, 2 * k - 1, 1. / k / ( k - 1 )});

      // Singlet
      _rules[1] = { {M4, 0, 1}, {M0, 1, 1}, {M5, 1, 1./6.}, {M6, onf, - 1. / nf1} };
      for (int k = nf + 2; k <= 6; k++)
        _rules[1].push_back({M5, 2 * k - 1, 1. / k / ( k - 1 )});

      // Light singlet-like distributions
      for (int k = 2; k <= nf; k++)
        _rules[2 * k - 1] = { {M0, 2 * k - 1, 1}, {M7, 2 * k - 1, 1} };

      // Heavy singlet-like distributions
      _rules[onf] = { {M4, 0, - static_cast<double>(nf)}, {M5, 1, - nf / 6.}, {M7, 1, nf * nf1 / 6.}, {M0, onf, 1}, {M6, onf, static_cast<double>(nf) / nf1}, {M7, onf, 1} };
      for (int k = nf + 2; k <= 6; k++)
        {
          _rules[onf].push_back({M5, 2 * k - 1, - nf / (static_cast<double>(k * ( k - 1 )))});
          _rules[onf].push_back({M7, 2 * k - 1,  nf * nf1 / (static_cast<double>(k * ( k - 1 )))});
        }

      // Super-heavy singlet-like distributions. They match exactly
      // like the singlet.
      for (int l = nf + 2; l <= 6; l++)
        {
          _rules[2 * l - 1] = { {M4, 0, 1}, {M5, 1, 1./6.}, {M6, onf, - 1. / nf1} };
          for (int k = nf + 2; k <= 6; k++)
            {
              if (k == l)
                _rules[2 * l - 1].push_back({M0, 2 * k - 1, 1});
              _rules[2 * l - 1].push_back({M5, 2 * k - 1, 1. / k / ( k - 1 )});
            }
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
   * ConvolutionMap specialised for the matching of operators using
   * the QCD evolution basis.
   */
  class MatchingOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PNS, PQQ, PQG, PGQ, PGG, PT3Q, PT8Q, PT15Q, PT24Q, PT35Q, PT3G, PT8G, PT15G, PT24G, PT35G};
    enum Object:  int {GG, GQ, QG, QQ, VAL, T3S, T3G, V3V, T8S, T8G, V8V, T15S, T15G, V15V, T24S, T24G, V24V, T35S, T35G, V35V};

    /**
     * @brief The MatchingOperatorBasisQCD constructor for the
     * matching in the QCD evolution basis with nf active flavours.
     * @param nf: number of active flavours
     */
    MatchingOperatorBasisQCD(int const& nf):
      ConvolutionMap{"MatchingOperatorBasisQCD_" + std::to_string(nf)}
    {
      _rules[GG]  = { {PGG, GG, 1}, {PGQ, QG, 1} };
      _rules[GQ]  = { {PGG, GQ, 1}, {PGQ, QQ, 1} };
      _rules[QG]  = { {PQG, GG, 1}, {PQQ, QG, 1} };
      _rules[QQ]  = { {PQG, GQ, 1}, {PQQ, QQ, 1} };
      _rules[VAL] = { {PNS, VAL, 1} };
      for (int k = 1; k < 6; k++)
        if (k < nf)
          {
            _rules[3 * k + 2] = { {PNS, 3 * k + 2, 1} };
            _rules[3 * k + 3] = { {PNS, 3 * k + 3, 1} };
            _rules[3 * k + 4] = { {PNS, 3 * k + 4, 1} };
          }
        else
          {
            _rules[3 * k + 2] = { {4 + k, QQ, 1}, {9 + k, GQ, 1} };
            _rules[3 * k + 3] = { {4 + k, QG, 1}, {9 + k, GG, 1} };
            _rules[3 * k + 4] = _rules[VAL];
          }
    };
  };
  ///@}
}
