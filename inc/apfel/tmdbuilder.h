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
    double                    Threshold;
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
   * @name TMD object initializers
   * Collection of functions that initialise TmdObjects structure for
   * the perturbartive evolution and matching currently available.
   */
  ///@{
  /**
   * @brief The InitializeDglapObjectsQCD function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   */

  /**
   * @brief The ComputeMatchingFunctions precomputes the constant
   * (i.e. non-logaritmic) perturbative coefficients of the matching
   * functions to match PDFs/FFs on the respective TMDs at small
   * values of bT as well as the perturbative coefficients of the
   * anomalous dimensions necessary for the Collins-Soper evolution.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map of TmdObjects objects, one for each possible nf
   */
    map<int,TmdObjects> InitializeTmdObjects(Grid           const& g,
					     vector<double> const& Thresholds,
					     double         const& IntEps = 1e-5);
  ///@}

  /**
   * @name TMD builders
   * Collection of functions that build a TMD distributions (both PDFs
   * and FF) as Set<Distribution>-valued functions. These functions
   * perform evolution and matching either separately or alltogether.
   */
  ///@{
  /**
   * @brief Function that returns the matched and evolved TMD PDFs in
   * b-space as functions of the final scale and rapidity. It assumes
   * the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param DglapObj: the (space-like) DGLAP objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param fNP: the non-perturbative function
   * @param Mu0b: the matching scale function
   * @param Mub: the initial scale function
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved TMD PDFs
   */
    function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(map<int,TmdObjects>                            const& TmdObj,
											  map<int,DglapObjects>                          const& DglapObj,
											  function<Set<Distribution>(double const&)>     const& CollPDFs,
											  function<double(double const&, double const&)> const& fNP,
											  function<double(double const&)>                const& Mu0b,
											  function<double(double const&)>                const& Mub,
											  int                                            const& PerturbativeOrder,
											  function<double(double const&)>                const& Alphas,
											  double                                         const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched and evolved TMD FFs in
   * b-space as functions of the final scale and rapidity. It assumes
   * the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param DglapObj: the (time-like) DGLAP objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param fNP: the non-perturbative function
   * @param Mu0b: the matching scale function
   * @param Mub: the initial scale function
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved TMD FFs
   */
    function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdFFs(map<int,TmdObjects>                            const& TmdObj,
											 map<int,DglapObjects>                          const& DglapObj,
											 function<Set<Distribution>(double const&)>     const& CollFFs,
											 function<double(double const&, double const&)> const& fNP,
											 function<double(double const&)>                const& Mu0b,
											 function<double(double const&)>                const& Mub,
											 int                                            const& PerturbativeOrder,
											 function<double(double const&)>                const& Alphas,
											 double                                         const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched TMD PDFs in b-space. It
   * assumes the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param DglapObj: the (space-like) DGLAP objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param fNP: the non-perturbative function
   * @param Mub: the initial scale function
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD PDFs
   */
    function<Set<Distribution>(double const&)> MatchTmdPDFs(map<int,TmdObjects>                        const& TmdObj,
							    map<int,DglapObjects>                      const& DglapObj,
							    function<Set<Distribution>(double const&)> const& CollPDFs,
							    function<double(double const&)>            const& Mub,
							    int                                        const& PerturbativeOrder,
							    function<double(double const&)>            const& Alphas);

  /**
   * @brief Function that returns the matched TMD FFs in b-space. It
   * assumes the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param DglapObj: the (time-like) DGLAP objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param fNP: the non-perturbative function
   * @param Mub: the initial scale function
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD FFs
   */
    function<Set<Distribution>(double const&)> MatchTmdFFs(map<int,TmdObjects>                        const& TmdObj,
							   map<int,DglapObjects>                      const& DglapObj,
							   function<Set<Distribution>(double const&)> const& CollPDFs,
							   function<double(double const&)>            const& Mub,
							   int                                        const& PerturbativeOrder,
							   function<double(double const&)>            const& Alphas);

  /**
   * @brief Function that returns the evolution factors for gluon and quarks. It
   * assumes the zeta-prescription. Universal for PDFs and FFs.
   * @param TmdObj: the TMD objects
   * @param Mu0b: the matching scale function
   * @param Mub: the initial scale function
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return vector<double>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. The 0-th component contains the
   * gluon evolution vactor, the remaining 12, from 1 to 12, are all
   * equal and represent the quark evolution factors.
   */
  //_____________________________________________________________________________
  function<vector<double>(double const&, double const&, double const&)> EvolutionFactors(map<int,TmdObjects>             const& TmdObj,
											 function<double(double const&)> const& Mu0b,
											 function<double(double const&)> const& Mub,
											 int                             const& PerturbativeOrder,
											 function<double(double const&)> const& Alphas,
											 double                          const& IntEps = 1e-7);
  ///@}
}
