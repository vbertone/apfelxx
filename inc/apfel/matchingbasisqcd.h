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
    enum Operand: int {M0, M1, M2, M3, M4, M5, M6, M7, M8, M9};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The MatchingBasisQCD constructor for the matching in the
     * QCD evolution basis with nf active flavours.
     @param nf: number of active flavours
     */
    MatchingBasisQCD(int const& nf);
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
    MatchingOperatorBasisQCD(int const& nf);
  };
  ///@}
}
