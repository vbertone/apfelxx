//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/lagrangeinterpolator.h>

namespace apfel
{
  /**
   * @brief The Distribution class for PDFs.
   *
   * This class provides methods to inherit a custom PDF distribution
   * in a generic basis.
   */
  class Distribution: public LagrangeInterpolator
  {
  public:
    /**
     * @brief Distribution constructors.
     */
    Distribution(Grid const& gr);
  };

}
