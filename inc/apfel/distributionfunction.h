//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/distribution.h"

#include <functional>

using std::function;

namespace apfel {

  /**
   * @brief The PDF class
   *
   * Helper class for the construction of PDFs from function
   */
  class DistributionFunction: public Distribution
  {
  public:

    /**
     * @brief Default constructor.
     */
    DistributionFunction(Grid                                        const& g,
			 function<double(int const&, double const&)> const& InPDFsFunc,
			 int                                         const& ipdf);

    /**
     * @brief Default constructor.
     */
    DistributionFunction(Grid                                                       const& g,
			 function<double(int const&, double const&, double const&)> const& InPDFsFunc,
			 int                                                        const& ipdf,
			 double                                                     const& Q);

  };

}
