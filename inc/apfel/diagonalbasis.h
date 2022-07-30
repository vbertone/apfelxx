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
   * @defgroup DiagonalBasis Diagonal convolution maps
   * Simple diagonal convolution basis
   * @ingroup ConvMap
   */
  ///@{
  /**
   * @brief The DiagonalBasis class is the simplest derivation of
   * ConvolutionMap meant to essentially perform a scalar product of
   * two sets of objects.
   */
  class DiagonalBasis: public ConvolutionMap
  {
  public:
    /**
     * @brief The DiagonalBasis constructor
     * @param nf: number of elements
     * @param offset: starting index for the enumeration on the distributions (default: 0)
     */
    DiagonalBasis(int const& nf, int const& offset = 0);
  };
  ///@}
}
