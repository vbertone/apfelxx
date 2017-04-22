//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/observable.h"

#include <functional>
#include <vector>

using std::function;
using std::vector;

namespace apfel
{
  /**
   * @brief The F2BuildZM
   */
  unordered_map<int,Observable> F2BuildZM(Grid                                                       const& g,
					  function<double(int const&, double const&, double const&)> const& InDistFunc,
					  vector<double>                                             const& Thresholds,
					  int                                                        const& PerturbativeOrder,
					  function<double(double const&)>                            const& Alphas,
					  function<vector<double>(double const&)>                    const& Charges,
					  double                                                     const& IntEps = 1e-5);

  /**
   * @brief The FLBuildZM
   */
  unordered_map<int,Observable> FLBuildZM(Grid                                                       const& g,
					  function<double(int const&, double const&, double const&)> const& InDistFunc,
					  vector<double>                                             const& Thresholds,
					  int                                                        const& PerturbativeOrder,
					  function<double(double const&)>                            const& Alphas,
					  function<vector<double>(double const&)>                    const& Charges,
					  double                                                     const& IntEps = 1e-5);

  /**
   * @brief The F3BuildZM
   */
  unordered_map<int,Observable> F3BuildZM(Grid                                                       const& g,
					  function<double(int const&, double const&, double const&)> const& InDistFunc,
					  vector<double>                                             const& Thresholds,
					  int                                                        const& PerturbativeOrder,
					  function<double(double const&)>                            const& Alphas,
					  function<vector<double>(double const&)>                    const& Charges,
					  double                                                     const& IntEps = 1e-5);
}
