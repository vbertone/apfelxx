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
    map<int,Set<Operator>> C0;
    map<int,Set<Operator>> C1;
    map<int,map<int,Set<Operator>>> C2;
    map<int,ConvolutionMap> ConvBasis;
    function<ConvolutionMap(vector<double> const&)> ConvBasisTot;
  };

  /**
   * @brief The InitializeF2NCObjectsMassive, precompute the
   * perturbative coefficients of coefficient functions as functions
   * of eta(Q) for the massive NC structure function F2 and store them
   * in a function of eta that returns a 'StructureFunctionObjects'
   * structure.
   *
   * @param g the grid
   * @param IntEps the integration accuracy
   * @return
   */
  function<StructureFunctionObjects(double const&)> InitializeF2NCObjectsMassive(Grid           const& g,
										 vector<double> const& Masses,
										 double         const& IntEps = 1e-5,
										 int            const& nxi    = 150,
										 double         const& ximin  = 0.001,
										 double         const& ximax  = 100000,
										 int            const& intdeg = 3,
										 double         const& lambda = 0.0005);

  /**
   * @brief The InitializeF2ObjectsZM, precompute the perturbative
   * coefficients of coefficient functions for NC F2 and store them in
   * the 'StructureFunctionObjects' structure.
   *
   * @param g the grid
   * @param IntEps the integration accuracy
   * @return
   */
  StructureFunctionObjects InitializeF2NCObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for FL NC.
   */
  StructureFunctionObjects InitializeFLNCObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for F3 NC.
   */
  StructureFunctionObjects InitializeF3NCObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for ( F2(nu) + F2(nubar) ) / 2 CC.
   */
  StructureFunctionObjects InitializeF2CCPlusObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for ( F2(nu) - F2(nubar) ) / 2 CC.
   */
  StructureFunctionObjects InitializeF2CCMinusObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for ( FL(nu) + FL(nubar) ) / 2 CC.
   */
  StructureFunctionObjects InitializeFLCCPlusObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for ( FL(nu) - FL(nubar) ) / 2 CC.
   */
  StructureFunctionObjects InitializeFLCCMinusObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for ( F3(nu) + F3(nubar) ) / 2 CC.
   */
  StructureFunctionObjects InitializeF3CCPlusObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief Same as above for ( F3(nu) - F3(nubar) ) / 2 CC.
   */
  StructureFunctionObjects InitializeF3CCMinusObjectsZM(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief The StructureFunctionBuildNC class constructs a map of
   * "Observable" objects.
   */
  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(StructureFunctionObjects                                const& FObj,
					      function<map<int,double>(double const&, double const&)> const& InDistFunc,
					      vector<double>                                          const& Thresholds,
					      int                                                     const& PerturbativeOrder,
					      function<double(double const&)>                         const& Alphas,
					      function<vector<double>(double const&)>                 const& Charges);

  /**
   * @brief Same as above but using a function of ipdf, x, and Q as an
   * input.
   */
  map<int,Observable> BuildStructureFunctions(StructureFunctionObjects                                   const& FObj,
					      function<double(int const&, double const&, double const&)> const& InDistFunc,
					      vector<double>                                             const& Thresholds,
					      int                                                        const& PerturbativeOrder,
					      function<double(double const&)>                            const& Alphas,
					      function<vector<double>(double const&)>                    const& Charges);

  /**
   * @brief The StructureFunctionBuildNC class constructs a map of
   * "Observable" objects. As compared to the functions above, this
   * takes as an input a function of Q as
   * StructureFunctionObjects. This is needed for th massive structure
   * functions.
   */
  //_____________________________________________________________________________
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&)>       const& FObj,
					      function<map<int,double>(double const&, double const&)> const& InDistFunc,
					      vector<double>                                          const& Thresholds,
					      int                                                     const& PerturbativeOrder,
					      function<double(double const&)>                         const& Alphas,
					      function<vector<double>(double const&)>                 const& Charges);

  /**
   * @brief Same as above but using a function of ipdf, x, and Q as an
   * input.
   */
  map<int,Observable> BuildStructureFunctions(function<StructureFunctionObjects(double const&)>          const& FObj,
					      function<double(int const&, double const&, double const&)> const& InDistFunc,
					      vector<double>                                             const& Thresholds,
					      int                                                        const& PerturbativeOrder,
					      function<double(double const&)>                            const& Alphas,
					      function<vector<double>(double const&)>                    const& Charges);

}
