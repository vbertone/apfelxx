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
    double                              Threshold;
    std::map<int,double>                Beta;
    std::map<int,double>                GammaCuspq;
    std::map<int,double>                GammaCuspg;
    std::map<int,double>                GammaVq;
    std::map<int,double>                GammaVg;
    std::map<int,std::valarray<double>> CSdq;
    std::map<int,std::valarray<double>> CSdg;
    std::map<int,std::valarray<double>> Lzetaq;
    std::map<int,std::valarray<double>> Lzetag;
    std::map<int,Set<Operator>>         MatchingFunctionsPDFs;
    std::map<int,Set<Operator>>         MatchingFunctionsFFs;
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
    std::map<int,TmdObjects> InitializeTmdObjects(Grid                const& g,
						  std::vector<double> const& Thresholds,
						  double              const& IntEps = 1e-5);
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
    std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(std::map<int,TmdObjects>                            const& TmdObj,
											       std::map<int,DglapObjects>                          const& DglapObj,
											       std::function<Set<Distribution>(double const&)>     const& CollPDFs,
											       std::function<double(double const&, double const&)> const& fNP,
											       std::function<double(double const&)>                const& Mu0b,
											       std::function<double(double const&)>                const& Mub,
											       int                                                 const& PerturbativeOrder,
											       std::function<double(double const&)>                const& Alphas,
											       double                                              const& IntEps = 1e-7);

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
    std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdFFs(std::map<int,TmdObjects>                            const& TmdObj,
											      std::map<int,DglapObjects>                          const& DglapObj,
											      std::function<Set<Distribution>(double const&)>     const& CollFFs,
											      std::function<double(double const&, double const&)> const& fNP,
											      std::function<double(double const&)>                const& Mu0b,
											      std::function<double(double const&)>                const& Mub,
											      int                                                 const& PerturbativeOrder,
											      std::function<double(double const&)>                const& Alphas,
											      double                                              const& IntEps = 1e-7);

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
    std::function<Set<Distribution>(double const&)> MatchTmdPDFs(std::map<int,TmdObjects>                        const& TmdObj,
								 std::map<int,DglapObjects>                      const& DglapObj,
								 std::function<Set<Distribution>(double const&)> const& CollPDFs,
								 std::function<double(double const&)>            const& Mub,
								 int                                             const& PerturbativeOrder,
								 std::function<double(double const&)>            const& Alphas);

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
    std::function<Set<Distribution>(double const&)> MatchTmdFFs(std::map<int,TmdObjects>                        const& TmdObj,
								std::map<int,DglapObjects>                      const& DglapObj,
								std::function<Set<Distribution>(double const&)> const& CollPDFs,
								std::function<double(double const&)>            const& Mub,
								int                                             const& PerturbativeOrder,
								std::function<double(double const&)>            const& Alphas);

  /**
   * @brief Function that returns the evolution factors for gluon and quarks. It
   * assumes the zeta-prescription. Universal for PDFs and FFs.
   * @param TmdObj: the TMD objects
   * @param Mu0b: the matching scale function
   * @param Mub: the initial scale function
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the strong coupling function
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return std::vector<double>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. The 0-th component contains the
   * gluon evolution vactor, the remaining 12, from 1 to 12, are all
   * equal and represent the quark evolution factors.
   */
  //_____________________________________________________________________________
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactors(std::map<int,TmdObjects>             const& TmdObj,
												   std::function<double(double const&)> const& Mu0b,
												   std::function<double(double const&)> const& Mub,
												   int                                  const& PerturbativeOrder,
												   std::function<double(double const&)> const& Alphas,
												   double                               const& IntEps = 1e-7);
  ///@}
}
