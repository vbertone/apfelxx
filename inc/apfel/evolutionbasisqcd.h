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
   * @brief The map between pair of indices corresponding to the
   * position of the operator in the evolution matrix and its linear
   * index.
   */
  // *INDENT-OFF*
  const std::map<std::pair<int, int>, int> Gkj =
    {
      //    g         Sigma        V           T3         V3          T8         V8         T15         V15        T24          V24          T35          V35
      {{ 0,0}, 0}, {{ 0,1}, 1},           {{ 0,3}, 2},           {{ 0,5}, 3},           {{ 0,7}, 4},           {{ 0,9}, 5},             {{ 0,11}, 6},             // g
      {{ 1,0}, 7}, {{ 1,1}, 8},           {{ 1,3}, 9},           {{ 1,5},10},           {{ 1,7},11},           {{ 1,9},12},             {{ 1,11},13},             // Sigma
                               {{2,2},14},                                                                                                                        // V
      {{ 3,0},15}, {{ 3,1},16},           {{ 3,3},17},           {{ 3,5},18},           {{ 3,7},19},           {{ 3,9},20},             {{ 3,11},21},             // T3
                                                      {{4,4},22},                                                                                                 // V3
      {{ 5,0},23}, {{ 5,1},24},           {{ 5,3},25},           {{ 5,5},26},           {{ 5,7},27},           {{ 5,9},28},             {{ 5,11},29},             // T8
                                                                             {{6,6},30},                                                                          // V8
      {{ 7,0},31}, {{ 7,1},32},           {{ 7,3},33},           {{ 7,5},34},           {{ 7,7},35},           {{ 7,9},36},             {{ 7,11},37},             // T15
                                                                                                    {{8,8},38},                                                   // V15
      {{ 9,0},39}, {{ 9,1},40},           {{ 9,3},41},           {{ 9,5},42},           {{ 9,7},43},           {{ 9,9},44},             {{ 9,11},45},             // T24
                                                                                                                           {{10,10},46},                          // V24
      {{11,0},47}, {{11,1},48},           {{11,3},49},           {{11,5},50},           {{11,7},51},           {{11,9},52},             {{11,11},53},             // T35
                                                                                                                                                     {{12,12},54} // V35
    };
  // *INDENT-ON*

  /**
   * @defgroup ConvMap Convolution maps
   * Collection of convolution maps to combine sets of objects
   * according to the task.
   */
  ///@{
  ///@}

  /**
   * @defgroup EvolBases Evolution convolution maps
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for the DGLAP evolution in the VFNS.
   * @ingroup ConvMap
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
    EvolutionBasisQCD(int const& nf);
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
     * @brief The EvolutionOperatorBasisQCD constructor
     * @param nf: number of active flavours
     */
    EvolutionOperatorBasisQCD(int const& nf);
  };

  /**
   * @brief The EvolveDistributionsBasisQCD class is a derived of
   * ConvolutionMap specialised to match a set of evolution operators
   * to a set a initial-scale distributions.
   */
  class EvolveDistributionsBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The EvolveDistributionsBasisQCD constructor
     */
    EvolveDistributionsBasisQCD();
  };
  ///@}
}
