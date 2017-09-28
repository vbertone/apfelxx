//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/operator.h"
#include "apfel/set.h"

#include <valarray>

using std::valarray;

namespace apfel
{
  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to perform the TMD evolution and matchinf to the collinear
   * PDFs, i.e. perturbative coefficients of matching functions and
   * and all anomalous dimensions.
   */
  struct TmdObjects
  {
    map<int,double>           Beta;
    map<int,double>           GammaCuspq;
    map<int,double>           GammaCuspg;
    map<int,double>           GammaVq;
    map<int,double>           GammaVg;
    map<int,valarray<double>> CSdq;
    map<int,valarray<double>> CSdg;
    map<int,Set<Operator>>    MatchingFunctions;
  };

  /**
   * @brief The ComputeMatchingFunctions precomputes the constant
   * (i.e. non-logaritmic) perturbative coefficients of the matching
   * functions to match PDFs/FFs on the respective TMDs at small
   * values of bT.
   * @param g the grid
   * @param IntEps the integration accuracy
   * @return
   */
    map<int,TmdObjects> InitializeTmdObjects(Grid           const& g,
					     vector<double> const& Thresholds,			      
					     double         const& IntEps = 1e-5);
}
