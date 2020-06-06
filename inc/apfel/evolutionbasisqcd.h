//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

namespace apfel
{
  /**
   * @defgroup EvolBases Evolution convolution maps Collection of
   * derived classes from ConvolutionMap that implement the
   * convolution map for the DGLAP evolution in the VFNS.
   */
  ///@{
  /**
   * @brief The EvolutionBasisQCD class is a derived of ConvolutionMap
   * specialised for the DGLAP evolution of distributions using the
   * QCD evolution basis.
   */
  class EvolutionBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PNSP, PNSM, PNSV, PQQ, PQG, PGQ, PGG};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The EvolutionBasisQCD constructor for the DGLAP
     * evolution in the QCD evolution basis with nf active flavours.
     * @param nf: number of active flavours
     */
    EvolutionBasisQCD(int const& nf):
      ConvolutionMap{"EvolutionBasisQCD_" + std::to_string(nf)}
    {
      _rules[GLUON] = { {PGG, GLUON, 1}, {PGQ, SIGMA, static_cast<double>(nf) / 6} };
      for (int j = nf + 1; j <= 6; j++)
        _rules[GLUON].push_back({PGQ, 2 * j - 1, static_cast<double>(nf) / j / ( j - 1 )});

      _rules[SIGMA]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, static_cast<double>(nf) / 6} };
      for (int j = nf + 1; j <= 6; j++)
        _rules[SIGMA].push_back({PQQ, 2 * j - 1, static_cast<double>(nf) / j / ( j - 1 )});

      _rules[VALENCE] = { {PNSV, VALENCE, 1} };

      for (int i = 2; i <= nf; i++)
        {
          _rules[2 * i - 1] = { {PNSP, 2 * i - 1, 1} };
          _rules[2 * i]     = { {PNSM, 2 * i, 1} };
        }

      for (int i = nf + 1; i <= 6; i++)
        {
          _rules[2 * i - 1]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, static_cast<double>(nf) / 6} };
          for (int j = nf + 1; j <= 6; j++)
            _rules[2 * i - 1].push_back({PQQ, 2 * j - 1, static_cast<double>(nf) / j / ( j - 1 )});
          _rules[2 * i]     = { {PNSV, 2 * i, 1} };
        }
    };
  };

  /**
   * @brief The EvolutionOperatorBasisQCD class is a derived of
   * ConvolutionMap specialised for the DGLAP evolution of operators
   * using the QCD evolution basis.
   */
  class EvolutionOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PNSP, PNSM, PNSV, PQQ, PQG, PGQ, PGG};
    enum Object:  int {GG, GQ, QG, QQ, VAL, T3S, T3G, V3V, T8S, T8G, V8V, T15S, T15G, V15V, T24S, T24G, V24V, T35S, T35G, V35V};

    /**
     * @brief The EvolutionOperatorBasisQCD constructor for the DGLAP
     * evolution in the QCD evolution basis with nf active flavours.
     * @param nf: number of active flavours
     */
    EvolutionOperatorBasisQCD(int const& nf):
      ConvolutionMap{"EvolutionOperatorBasisQCD_" + std::to_string(nf)}
    {
      _rules[GG]  = { {PGG, GG, 1}, {PGQ, QG, 1} };
      _rules[GQ]  = { {PGG, GQ, 1}, {PGQ, QQ, 1} };
      _rules[QG]  = { {PQG, GG, 1}, {PQQ, QG, 1} };
      _rules[QQ]  = { {PQG, GQ, 1}, {PQQ, QQ, 1} };
      _rules[VAL] = { {PNSV, VAL, 1} };
      for (int k = 1; k < 6; k++)
        if (k < nf)
          {
            _rules[3 * k + 2] = { {PNSP, 3 * k + 2, 1} };
            _rules[3 * k + 3] = { {PNSP, 3 * k + 3, 1} };
            _rules[3 * k + 4] = { {PNSM, 3 * k + 4, 1} };
          }
        else
          {
            _rules[3 * k + 2] = { {PQG, GQ, 1}, {PQQ, 3 * k + 2, 1} };
            _rules[3 * k + 3] = { {PQG, GG, 1}, {PQQ, 3 * k + 3, 1} };
            _rules[3 * k + 4] = { {PNSV, 3 * k + 4, 1} };
          }
    };
  };

  /**
   * @brief The MatchEvolOperatorBasisQCD class is a derived of
   * ConvolutionMap specialised for the matching of operators at the
   * heavy-quark thresholds using the QCD evolution basis.
   */
  class MatchEvolOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {GG, GQ, QG, QQ, VAL, T3S, T3G, V3V, T8S, T8G, V8V, T15S, T15G, V15V, T24S, T24G, V24V, T35S, T35G, V35V};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The MatchEvolOperatorBasisQCD constructor for the
     * matching in the QCD evolution basis with nf active flavours.
     * @param nf: number of active flavours
     */
    MatchEvolOperatorBasisQCD(int const& nf):
      ConvolutionMap{"MatchEvolOperatorBasisQCD_" + std::to_string(nf)}
    {
      _rules[GLUON]   = { {GG, GLUON, 1}, {GQ, SIGMA, 1} };
      _rules[SIGMA]   = { {QG, GLUON, 1}, {QQ, SIGMA, 1} };
      _rules[VALENCE] = { {VAL, VALENCE, 1} };
      for (int k = 1; k < 6; k++)
        if (k < nf)
          {
            _rules[2 * k + 1] = { {3 * k + 2, 2 * k + 1, 1} };
            _rules[2 * k + 2] = { {3 * k + 4, 2 * k + 2, 1} };
          }
        else
          {
            _rules[2 * k + 1] = { {3 * k + 2, SIGMA, 1}, {3 * k + 3, GLUON, 1} };
            _rules[2 * k + 2] = { {3 * k + 4, VALENCE, 1} };
          }
    };
  };
  ///@}
}
