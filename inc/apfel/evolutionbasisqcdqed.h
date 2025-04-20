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
   * index in the evolution basis used for QCDxQED evolution.
   */
  // *INDENT-OFF*
  const std::map<std::pair<int, int>, int> GkjQCDQED =
    {
      //tau^+        mu^+           e^+            t^+            c^+            u^+            b^+            s^+            d^+            g              gamma
      {{0,0},0},     {{0,1},1},     {{0,2},2},     {{0,3},3},     {{0,4},4},     {{0,5},5},     {{0,6},6},     {{0,7},7},     {{0,8},8},     {{0,9},9},     {{0,10},10},   // tau^+
      {{1,0},11},    {{1,1},12},    {{1,2},13},    {{1,3},14},    {{1,4},15},    {{1,5},16},    {{1,6},17},    {{1,7},18},    {{1,8},19},    {{1,9},20},    {{1,10},21},   // mu^+
      {{2,0},22},    {{2,1},23},    {{2,2},24},    {{2,3},25},    {{2,4},26},    {{2,5},27},    {{2,6},28},    {{2,7},29},    {{2,8},30},    {{2,9},31},    {{2,10},32},   // e^+
      {{3,0},33},    {{3,1},34},    {{3,2},35},    {{3,3},36},    {{3,4},37},    {{3,5},38},    {{3,6},39},    {{3,7},40},    {{3,8},41},    {{3,9},42},    {{3,10},43},   // t^+
      {{4,0},44},    {{4,1},45},    {{4,2},46},    {{4,3},47},    {{4,4},48},    {{4,5},49},    {{4,6},50},    {{4,7},51},    {{4,8},52},    {{4,9},53},    {{4,10},54},   // c^+
      {{5,0},55},    {{5,1},56},    {{5,2},57},    {{5,3},58},    {{5,4},59},    {{5,5},60},    {{5,6},61},    {{5,7},62},    {{5,8},63},    {{5,9},64},    {{5,10},65},   // u^+
      {{6,0},66},    {{6,1},67},    {{6,2},68},    {{6,3},69},    {{6,4},70},    {{6,5},71},    {{6,6},72},    {{6,7},73},    {{6,8},74},    {{6,9},75},    {{6,10},76},   // b^+
      {{7,0},77},    {{7,1},78},    {{7,2},79},    {{7,3},80},    {{7,4},81},    {{7,5},82},    {{7,6},83},    {{7,7},84},    {{7,8},85},    {{7,9},86},    {{7,10},87},   // s^+
      {{8,0},88},    {{8,1},89},    {{8,2},90},    {{8,3},91},    {{8,4},92},    {{8,5},93},    {{8,6},94},    {{8,7},95},    {{8,8},96},    {{8,9},97},    {{8,10},98},   // d^+
      {{9,0},99},    {{9,1},100},   {{9,2},101},   {{9,3},102},   {{9,4},103},   {{9,5},104},   {{9,6},105},   {{9,7},106},   {{9,8},107},   {{9,9},108},   {{9,10},109},  // g
      {{10,0},110},  {{10,1},111},  {{10,2},112},  {{10,3},113},  {{10,4},114},  {{10,5},115},  {{10,6},116},  {{10,7},117},  {{10,8},118},  {{10,9},119},  {{10,10},120}, // gamma
      //d^-           s^-           b^-            u^-            c^-            t^-
      {{11,11},121}, {{11,12},122}, {{11,13},123}, {{11,14},124}, {{11,15},125}, {{11,16},126}, // d^-
      {{12,11},127}, {{12,12},128}, {{12,13},129}, {{12,14},130}, {{12,15},131}, {{12,16},132}, // s^-
      {{13,11},133}, {{13,12},134}, {{13,13},135}, {{13,14},136}, {{13,15},137}, {{13,16},138}, // b^-
      {{14,11},139}, {{14,12},140}, {{14,13},141}, {{14,14},142}, {{14,15},143}, {{14,16},144}, // u^-
      {{15,11},145}, {{15,12},146}, {{15,13},147}, {{15,14},148}, {{15,15},149}, {{15,16},150}, // c^-
      {{16,11},151}, {{16,12},152}, {{16,13},153}, {{16,14},154}, {{16,15},155}, {{16,16},156}, // t^-
      // e^-         mu^-           tau^-
      {{17,17},157},                               // e^-
                     {{18,18},158},                // mu^-
                                    {{19,19},159}  // tau^-
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
   * @brief The EvolutionBasisQCDQED class is a derived of ConvolutionMap
   * specialised for the DGLAP evolution of distributions using the
   * physical basis.
   */
  class EvolutionBasisQCDQED: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PPDD, PPUU, PPLL, PMDD, PMUU, PMLL,
                       PPSDD, PPSDU, PPSDL, PPSUD, PPSUU, PPSUL, PPSLD, PPSLU, PPSLL, PPV,
                       PDg, PUg, PDgm, PUgm, PLgm,
                       PgD, PgU, PgmD, PgmU, PgmL,
                       Pgg, Pggm, Pgmg, Pgmgm
                      };
    enum Object:  int {TAUP, MUP, EP, TP, CP, UP, BP, SP, DP, GLUON, PHOTON, DM, SM, BM, UM, CM, TM, EM, MUM, TAUM};

    /**
     * @brief The EvolutionBasisQCDQED constructor for the DGLAP
     * evolution in the evolution basis for QCDxQED evolution.
     * @param nd: number of active down-type quarks below threshold
     * @param nu: number of active up-type quarks below threshold
     * @param nl: number of active leptons below threshold
     */
    EvolutionBasisQCDQED(int const& nd, int const& nu, int const& nl);
  };
  ///@}
}
