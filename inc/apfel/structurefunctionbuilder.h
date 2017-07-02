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
#include "apfel/disbasis.h"

#include <functional>
#include <vector>

using std::function;
using std::vector;

namespace apfel
{

  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to compute the DIS structure functions, i.e. perturbative
   * coefficients of the coefficient functions.
   */
  struct StructureFunctionObjects
  {
    vector<int> skip;
    unordered_map<int,DISNCBasis> ConvBasis;
    unordered_map<int,Set<Operator>> C0;
    unordered_map<int,Set<Operator>> C1;
    unordered_map<int,unordered_map<int,Set<Operator>>> C2;
  };

  /**
   * @brief The InitializeF2ObjectsZM, precompute the perturbative
   * coefficients of coefficient functions for F2 and store them in
   * the 'StructureFunctionObjects' structure.
   *
   * @param g the grid
   * @param IntEps the integration accuracy
   * @return
   */
  StructureFunctionObjects InitializeF2ObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for FL.
   */
  StructureFunctionObjects InitializeFLObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for F3.
   */
  StructureFunctionObjects InitializeF3ObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief The F2BuildZM
   */
  //_____________________________________________________________________________
  unordered_map<int,Observable> StructureFunctionBuildNC(StructureFunctionObjects                                          const& FObj,
							 function<unordered_map<int,double>(double const&, double const&)> const& InDistFunc,
							 vector<double>                                                    const& Thresholds,
							 int                                                               const& PerturbativeOrder,
							 function<double(double const&)>                                   const& Alphas,
							 function<vector<double>(double const&)>                           const& Charges);

  /**
   * @brief The F2BuildZM
   */
  unordered_map<int,Observable> StructureFunctionBuildNC(StructureFunctionObjects                                   const& FObj,
							 function<double(int const&, double const&, double const&)> const& InDistFunc,
							 vector<double>                                             const& Thresholds,
							 int                                                        const& PerturbativeOrder,
							 function<double(double const&)>                            const& Alphas,
							 function<vector<double>(double const&)>                    const& Charges);

}
