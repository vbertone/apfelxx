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
   * @defgroup MatchBases Matching convolution maps
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for the matching at the heavy-quark
   * thresholds.
   */
  ///@{
  /**
   * @brief The MatchingBasisQCD class is a derived of ConvolutionMap
   * specialised for the matching of distributions using the QCD
   * evolution basis.
   */
  class MatchingBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PNS, PQQ, PQG, PGQ, PGG, PT3Q, PT8Q, PT15Q, PT24Q, PT35Q, PT3G, PT8G, PT15G, PT24G, PT35G};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The MatchingBasisQCD default constructor for the
     * matching in the QCD evolution basis with nf active flavours.
     * @param nf: number of active flavours
     */
  MatchingBasisQCD(int const& nf):
    ConvolutionMap{"MatchingBasisQCD_" + std::to_string(nf)}
    {
      _rules[GLUON] = { {PGG, GLUON, 1}, {PGQ, SIGMA, 1} };
      _rules[SIGMA] = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
      _rules[VALENCE] = { {PNS, VALENCE, 1} };
      for (int k = 1; k < 6; k++)
	if (k < nf)
	  {
	    _rules[2*k+1] = { {PNS, 2*k+1, 1} };
	    _rules[2*k+2] = { {PNS, 2*k+2, 1} };
	  }
	else
	  {
	    _rules[2*k+1] = { {4+k, SIGMA, 1}, {9+k, GLUON, 1} };
	    _rules[2*k+2] = _rules[VALENCE];
	  }
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
     * @brief The MatchingOperatorBasisQCD default constructor for the
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
	    _rules[3*k+2] = { {PNS, 3*k+2, 1} };
	    _rules[3*k+3] = { {PNS, 3*k+3, 1} };
	    _rules[3*k+4] = { {PNS, 3*k+4, 1} };
	  }
	else
	  {
	    _rules[3*k+2] = { {4+k, QQ, 1}, {9+k, GQ, 1} };
	    _rules[3*k+3] = { {4+k, QG, 1}, {9+k, GG, 1} };
	    _rules[3*k+4] = _rules[VAL];
	  }
    };
  };
  ///@}
}
