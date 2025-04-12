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
   * index in the physical basis.
   */
  // *INDENT-OFF*
  const std::map<std::pair<int, int>, int> GkjPhys =
    {
      //     t-                b-                c-                s-                u-                u-                 g                 d+                u+                s+                c+                b+                t+
      {{-6+6,-6+6}, 0}, {{-6+6,-5+6}, 1}, {{-6+6,-4+6}, 2}, {{-6+6,-3+6}, 3}, {{-6+6,-2+6}, 4}, {{-6+6,-1+6}, 5},                                                                                                                               // t-
      {{-5+6,-6+6}, 6}, {{-5+6,-5+6}, 7}, {{-5+6,-4+6}, 8}, {{-5+6,-3+6}, 9}, {{-5+6,-2+6},10}, {{-5+6,-1+6},11},                                                                                                                               // b-
      {{-4+6,-6+6},12}, {{-4+6,-5+6},13}, {{-4+6,-4+6},14}, {{-4+6,-3+6},15}, {{-4+6,-2+6},16}, {{-4+6,-1+6},17},                                                                                                                               // c-
      {{-3+6,-6+6},18}, {{-3+6,-5+6},19}, {{-3+6,-4+6},20}, {{-3+6,-3+6},21}, {{-3+6,-2+6},22}, {{-3+6,-1+6},23},                                                                                                                               // s-
      {{-2+6,-6+6},24}, {{-2+6,-5+6},25}, {{-2+6,-4+6},26}, {{-2+6,-3+6},27}, {{-2+6,-2+6},28}, {{-2+6,-1+6},29},                                                                                                                               // u-
      {{-1+6,-6+6},30}, {{-1+6,-5+6},31}, {{-1+6,-4+6},32}, {{-1+6,-3+6},33}, {{-1+6,-2+6},34}, {{-1+6,-1+6},35},                                                                                                                               // d-
                                                                                                                  {{ 0+6, 0+6},36}, {{ 0+6, 1+6},37}, {{ 0+6, 2+6},38}, {{ 0+6, 3+6},39}, {{ 0+6, 4+6},40}, {{ 0+6, 5+6},41}, {{ 0+6, 6+6},42}, // g
                                                                                                                  {{ 1+6, 0+6},43}, {{ 1+6, 1+6},44}, {{ 1+6, 2+6},45}, {{ 1+6, 3+6},46}, {{ 1+6, 4+6},47}, {{ 1+6, 5+6},48}, {{ 1+6, 6+6},49}, // d+
                                                                                                                  {{ 2+6, 0+6},50}, {{ 2+6, 1+6},51}, {{ 2+6, 2+6},52}, {{ 2+6, 3+6},53}, {{ 2+6, 4+6},54}, {{ 2+6, 5+6},55}, {{ 2+6, 6+6},56}, // u+
                                                                                                                  {{ 3+6, 0+6},57}, {{ 3+6, 1+6},58}, {{ 3+6, 2+6},59}, {{ 3+6, 3+6},60}, {{ 3+6, 4+6},61}, {{ 3+6, 5+6},62}, {{ 3+6, 6+6},63}, // s+
                                                                                                                  {{ 4+6, 0+6},64}, {{ 4+6, 1+6},65}, {{ 4+6, 2+6},66}, {{ 4+6, 3+6},67}, {{ 4+6, 4+6},68}, {{ 4+6, 5+6},69}, {{ 4+6, 6+6},70}, // c+
                                                                                                                  {{ 5+6, 0+6},71}, {{ 5+6, 1+6},72}, {{ 5+6, 2+6},73}, {{ 5+6, 3+6},74}, {{ 5+6, 4+6},75}, {{ 5+6, 5+6},76}, {{ 5+6, 6+6},77}, // b+
                                                                                                                  {{ 6+6, 0+6},78}, {{ 6+6, 1+6},79}, {{ 6+6, 2+6},80}, {{ 6+6, 3+6},81}, {{ 6+6, 4+6},82}, {{ 6+6, 5+6},83}, {{ 6+6, 6+6},84}  // y+
    };
  // *INDENT-ON*

  /**
   * @defgroup PhysBases Physical convolution maps
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for the DGLAP evolution in the VFNS in the
   * physical basis.
   * @ingroup ConvMap
   */
  ///@{
  /**
   * @brief The PhysicalBasisQCD class is a derived of ConvolutionMap
   * specialised for the DGLAP evolution of distributions using the
   * physical basis.
   */
  class PhysicalBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PNV, PPV, PNS, PPS, PQG, PGQ, PGG};
    enum Object:  int {TM, BM, CM, SM, UM, DM, GLUON, DP, UP, SP, CP, BP, TP};

    /**
     * @brief The PhysicalBasisQCD constructor for the DGLAP
     * evolution in the physical basis with nf active flavours.
     * @param nf: number of active flavours
     */
    PhysicalBasisQCD(int const& nf);
  };

  /**
   * @brief The PhysicalOperatorBasisQCD class is a derived of
   * ConvolutionMap specialised for the DGLAP evolution of operators
   * using the physical basis.
   */
  class PhysicalOperatorBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The PhysicalOperatorBasisQCD constructor
     * @param nf: number of active flavours
     */
    PhysicalOperatorBasisQCD(int const& nf);
  };

  /**
   * @brief The PhysicalEvolveDistributionsBasisQCD class is a derived
   * of ConvolutionMap specialised to match a set of evolution
   * operators to a set a initial-scale distributions in the physical
   * basis.
   */
  class PhysicalEvolveDistributionsBasisQCD: public ConvolutionMap
  {
  public:
    /**
     * @brief The PhysicalEvolveDistributionsBasisQCD constructor
     */
    PhysicalEvolveDistributionsBasisQCD();
  };
  ///@}
}
