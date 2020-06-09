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
      {{ 0,0}, 0}, {{ 0,1}, 1},           {{ 0,3}, 2},           {{ 0,5}, 3},           {{ 0,7}, 4},           {{ 0,9}, 5},             {{ 0,11}, 6},
      {{ 1,0}, 7}, {{ 1,1}, 8},           {{ 1,3}, 9},           {{ 1,5},10},           {{ 1,7},11},           {{ 1,9},12},             {{ 1,11},13},
                               {{2,2},14},
      {{ 3,0},15}, {{ 3,1},16},           {{ 3,3},17},           {{ 3,5},18},           {{ 3,7},19},           {{ 3,9},20},             {{ 3,11},21},
                                                      {{4,4},22},
      {{ 5,0},23}, {{ 5,1},24},           {{ 5,3},25},           {{ 5,5},26},           {{ 5,7},27},           {{ 5,9},28},             {{ 5,11},29},
                                                                             {{6,6},30},
      {{ 7,0},31}, {{ 7,1},32},           {{ 7,3},33},           {{ 7,5},34},           {{ 7,7},35},           {{ 7,9},36},             {{ 7,11},37},
                                                                                                    {{8,8},38},
      {{ 9,0},39}, {{ 9,1},40},           {{ 9,3},41},           {{ 9,5},42},           {{ 9,7},43},           {{ 9,9},44},             {{ 9,11},45},
                                                                                                                           {{10,10},46},
      {{11,0},47}, {{11,1},48},           {{11,3},49},           {{11,5},50},           {{11,7},51},           {{11,9},52},             {{11,11},53},
                                                                                                                                                     {{12,12},54}
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
    EvolutionBasisQCD(int const& nf):
      ConvolutionMap{"EvolutionBasisQCD_" + std::to_string(nf)}
    {
      _rules[GLUON] = { {PGG, GLUON, 1}, {PGQ, SIGMA, 1} };
      for (int j = nf + 1; j <= 6; j++)
        _rules[GLUON].push_back({PGQ, 2 * j - 1, 6. / j / ( j - 1 )});

      _rules[SIGMA]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
      for (int j = nf + 1; j <= 6; j++)
        _rules[SIGMA].push_back({PQQ, 2 * j - 1, 6. / j / ( j - 1 )});

      _rules[VALENCE] = { {PNSV, VALENCE, 1} };

      for (int i = 2; i <= nf; i++)
        {
          _rules[2 * i - 1] = { {PNSP, 2 * i - 1, 1} };
          _rules[2 * i]     = { {PNSM, 2 * i,     1} };
        }

      for (int i = nf + 1; i <= 6; i++)
        {
          _rules[2 * i - 1]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
          for (int j = nf + 1; j <= 6; j++)
            _rules[2 * i - 1].push_back({PQQ, 2 * j - 1, 6. / j / ( j - 1 )});
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
     * @brief The EvolutionOperatorBasisQCD constructor
     * @param nf: number of active flavours
     */
    EvolutionOperatorBasisQCD(int const& nf):
      ConvolutionMap{"EvolutionOperatorBasisQCD_" + std::to_string(nf)}
    {
      // Allocate EvolutionBasisQCD object to retrieve the splitting
      // matrix rules
      const EvolutionBasisQCD eb{nf};

      // Get matrix of coefficients
      const matrix<std::vector<double>> rc = eb.GetRuleMatrix();

      // Get matrix of operator indices
      const matrix<std::vector<int>> ri = eb.GetRuleIndices();

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
    EvolveDistributionsBasisQCD():
      ConvolutionMap{"EvolveDistributionsBasisQCD"}
    {
      // Construct set of rules
      for (int k = 0; k < 13; k++)
        for (int j = 0; j < 13; j++)
          {
            if (Gkj.count({k, j}) == 0)
              continue;
            _rules[k].push_back({Gkj.at({k, j}), j, 1});
          }
    };
  };
  ///@}
}
