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
   * @defgroup EvolBasisPlusMinus Evolution convolution maps in the plus-minus basis
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for the DGLAP evolution in the VFNS.
   * @ingroup ConvMap
   */
  ///@{
  /**
   * @brief The EvolBasisPlusMinus class is derives from
   * ConvolutionMap and implements a basis in which plus (q+qbar) and
   * minus (q-qbar) combinations are fully coupled.
   */
  class EvolBasisPlusMinus: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enums
     */
    enum Operand: int {PGG, PGD, PGU, PGS, PGC, PGB, PGT, PDG, PUG, PSG, PCG, PBG, PTG, PDD, PUU, PSS, PCC, PBB, PTT, PPS, PMP};
    enum Object:  int {GLUON, DWP, DWM, UPP, UPM, STP, STM, CHP, CHM, BTP, BTM, TPP, TPM};

    /**
     * @brief The EvolBasisPlusMinus constructor
     */
    EvolBasisPlusMinus();
  };
  ///@}
}
