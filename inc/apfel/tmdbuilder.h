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
#include "apfel/dglapbuilder.h"
#include "apfel/tabulateobject.h"

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
    map<int,valarray<double>> Lzetaq;
    map<int,valarray<double>> Lzetag;
    map<int,Set<Operator>>    MatchingFunctionsPDFs;
    map<int,Set<Operator>>    MatchingFunctionsFFs;
  };

  /**
   * @brief The ComputeMatchingFunctions precomputes the constant
   * (i.e. non-logaritmic) perturbative coefficients of the matching
   * functions to match PDFs/FFs on the respective TMDs at small
   * values of bT.
   * @param g the x grid
   * @param Thresholds
   * @param IntEps the integration accuracy
   * @return
   */
    map<int,TmdObjects> InitializeTmdObjects(Grid           const& g,
					     vector<double> const& Thresholds,
					     double         const& IntEps = 1e-5);

  /**
   * @brief Function that returns the evolved TMD PDFs in b-space as
   * functions of the final scale and rapidity zeta. It assumes the
   * zeta-prescription.
   */
    function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(map<int,TmdObjects>                            const& TmdObj,
											  map<int,DglapObjects>                          const& DglapObj,
											  TabulateObject<Set<Distribution>>              const& CollPDFs,
											  function<double(double const&, double const&)> const& fNP,
											  function<double(double const&)>                const& Mu0b,
											  function<double(double const&)>                const& Mub,
											  int                                            const& PerturbativeOrder,
											  function<double(double const&)>                const& Alphas,
											  double                                         const& IntEps = 1e-7);
}
