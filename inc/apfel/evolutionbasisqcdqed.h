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
   * @addtogroup EvolBases
   */
  ///@{
  /**
   * @brief The EvolutionBasisQCDQED class is a derived of
   * ConvolutionMap specialised for the DGLAP evolution of
   * distributions using the mixed QCDxQED evolution basis.
   */
  class EvolutionBasisQCDQED: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {PNSP, PNSM, PNSV,           //!< Non-singlet pure QCD
                       PGG, PGQ,                   //!< Singlet pure QCD
                       PQG, PQQ,
                       PNSPQED, PNSMQED,           //!< Non-singlet QED (including O(alpha_s * alpha) corrections)
                       PGGQED,  PGGMQED,  PGQQED,  //!< Singlet QED (including O(alpha_s * alpha) corrections)
                       PGMGQED, PGMGMQED, PGMQQED,
                       PQGQED,  PQGMQED,  PQQQED,
                       PGMLQED, PQLQED,            //!< Lepton component of the singlet QED
                       PLGMQED, PLLQED
                      };
    enum Object:  int {GLUON, GAMMA, SIGMA, DSIGMA, SIGMAL,
                       VALENCE, DVALENCE,
                       VALENCEL,
                       T1U, V1U, T2U, V2U,
                       T1D, V1D, T2D, V2D,
                       T1L, V1L, T2L, V2L
                      };

    /**
     * @brief The EvolutionBasisQCDQED constructor for the DGLAP
     * evolution in the QCD evolution basis with nf active flavours.
     * @param nu: number of active up-type quark flavours
     * @param nd: number of active down-type quark flavours
     * @param nl: number of active charged-lepton flavours
     */
    EvolutionBasisQCDQED(int const& nu, int const& nd, int const& nl);
  };
  ///@}
}
