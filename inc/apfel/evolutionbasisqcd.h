//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

using namespace std;

namespace apfel
{
  /**
   * @brief ConvolutionMap derivation for the QCD evolution basis.
   *
   * This class, following the derivation procedure from
   * ConvolutionMap implements the Basis enumerator with custom tags
   * for the objects.
   */
  class EvolutionBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enums
     */
    enum Operand: int {PNSP, PNSM, PNSV, PQQ, PQG, PGQ, PGG};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The class constructor
     */
  EvolutionBasisQCD(int const& nf):
    ConvolutionMap{"EvolutionBasisQCD_" + std::to_string(nf)}
    {
      _rules[GLUON]   = { {PGG, GLUON, 1}, {PGQ, SIGMA, 1} };
      _rules[SIGMA]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
      _rules[VALENCE] = { {PNSV, VALENCE, 1} };
      for (int k = 1; k < 6; k++)
	if (k < nf)
	  {
	    _rules[2*k+1] = { {PNSP, 2*k+1, 1} };
	    _rules[2*k+2] = { {PNSM, 2*k+2, 1} };
	  }
	else
	  {
	    _rules[2*k+1] = _rules[SIGMA];
	    _rules[2*k+2] = _rules[VALENCE];
	  }
    };
  };

  class EvolutionOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enums
     */
    enum Operand: int {PNSP, PNSM, PNSV, PQQ, PQG, PGQ, PGG};
    enum Object:  int {GG, GQ, QG, QQ, VAL, T3S, T3G, V3V, T8S, T8G, V8V, T15S, T15G, V15V, T24S, T24G, V24V, T35S, T35G, V35V};

    /**
     * @brief The class constructor
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
	    _rules[3*k+2] = { {PNSP, 3*k+2, 1} };
	    _rules[3*k+3] = { {PNSP, 3*k+3, 1} };
	    _rules[3*k+4] = { {PNSM, 3*k+4, 1} };
	  }
	else
	  {
	    _rules[3*k+2] = _rules[QQ];
	    _rules[3*k+3] = _rules[QG];
	    _rules[3*k+4] = _rules[VAL];
	  }
    };
  };

  class MatchEvolOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enums
     */
    enum Operand: int {GG, GQ, QG, QQ, VAL, T3S, T3G, V3V, T8S, T8G, V8V, T15S, T15G, V15V, T24S, T24G, V24V, T35S, T35G, V35V};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The class constructor
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
	    _rules[2*k+1] = { {3*k+2, 2*k+1, 1} };
	    _rules[2*k+2] = { {3*k+4, 2*k+2, 1} };
	  }
	else
	  {
	    _rules[2*k+1] = { {3*k+2, SIGMA, 1}, {3*k+3, GLUON, 1} };
	    _rules[2*k+2] = { {3*k+4, VALENCE, 1} };
	  }
    };
  };
}
