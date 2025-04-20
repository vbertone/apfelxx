//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"
#include "apfel/constants.h"

namespace apfel
{
  /**
   * @defgroup MatchBasesQCDQED Matching convolution maps for QCDxQED
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for the matching at the quark/lepton
   * thresholds.
   * @ingroup ConvMap
   */
  ///@{
  /**
   * @brief The MatchingBasisQCDQED class is a derived of
   * ConvolutionMap specialised for the matching of distributions in
   * QCDxQED evolutions.
   */
  class MatchingBasisQCDQED: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {ONE, KQg, KXgm, KQqp, KXX, Kqg, KNSq, Kqqp, Kgg, Kgq, KgQ, KXgmgm, KgmX};
    enum Object:  int {TAUP, MUP, EP, TP, CP, UP, BP, SP, DP, GLUON, PHOTON, DM, SM, BM, UM, CM, TM, EM, MUM, TAUM};

    /**
     * @brief The MatchingBasisQCDQED constructor for the
     * matching in the QCD evolution basis with nf active flavours.
     * @param nd: number of active down-type quarks below threshold
     * @param nu: number of active up-type quarks below threshold
     * @param nl: number of active leptons below threshold
     * @param species: parton species of the heavy parton
     */
    MatchingBasisQCDQED(int const& nd, int const& nu, int const& nl, PartonSpecies const& species);
  };

  /**
   * @brief The MatchingOperatorBasisQCDQED class is a derived of
   * ConvolutionMap specialised for the matching of the evolution of
   * operators at the heavy-quark thresholds in QCDxQED evolution.
   */
  class MatchingOperatorBasisQCDQED: public ConvolutionMap
  {
  public:
    /**
     * @brief The MatchingOperatorBasisQCDQED constructor
     * @param nd: number of active down-type quarks below threshold
     * @param nu: number of active up-type quarks below threshold
     * @param nl: number of active leptons below threshold
     * @param species: parton species of the heavy parton
     */
    MatchingOperatorBasisQCDQED(int const& nd, int const& nu, int const& nl, PartonSpecies const& species);
  };
  ///@}
}
