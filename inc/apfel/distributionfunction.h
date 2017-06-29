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
#include <unordered_map>

using std::function;
using std::unordered_map;

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
			 function<double(int const&, double const&)> const& InDistFunc,
			 int                                         const& ipdf);

    /**
     * @brief Default constructor.
     */
    DistributionFunction(Grid                                                       const& g,
			 function<double(int const&, double const&, double const&)> const& InDistFunc,
			 int                                                        const& ipdf,
			 double                                                     const& Q);


    /**
     * @brief Default constructor.
     */
    DistributionFunction(Grid                   const& g,
			 vector<double>         const& DistJointGrid,
			 vector<vector<double>> const& DistSubGrid);

  };

}
