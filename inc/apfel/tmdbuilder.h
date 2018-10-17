//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/dglapbuilder.h"
#include "apfel/tabulateobject.h"

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
    double                                   Threshold;
    std::map<int,double>                     Beta;
    std::map<int,double>                     GammaFq;
    std::map<int,double>                     GammaFg;
    std::map<int,double>                     GammaCusp;
    std::map<int,std::vector<double>>        GammaCS;
    std::map<int,std::vector<Set<Operator>>> MatchingFunctionsPDFs;
    std::map<int,std::vector<Set<Operator>>> MatchingFunctionsFFs;
  };

  /**
   * @name TMD object initializers
   * Collection of functions that initialise TmdObjects structure for
   * the perturbartive evolution and matching currently available.
   */
  ///@{
  /**
   * @brief The InitializeTmdObjects function precomputes the
   * perturbative coefficients required for the evolution and matching
   * of TMD PDFs and FFs and store them into a 'TmdObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
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
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved TMD PDFs
   */
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(std::map<int,TmdObjects>                        const& TmdObj,
											     std::function<Set<Distribution>(double const&)> const& CollPDFs,
											     std::function<double(double const&)>            const& Alphas,
											     int                                             const& PerturbativeOrder,
											     double                                          const& Ci = 1,
											     double                                          const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched and evolved TMD FFs in
   * b-space as functions of the final scale and rapidity. It assumes
   * the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved TMD FFs
   */
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdFFs(std::map<int,TmdObjects>                        const& TmdObj,
											    std::function<Set<Distribution>(double const&)> const& CollFFs,
											    std::function<double(double const&)>            const& Alphas,
											    int                                             const& PerturbativeOrder,
											    double                                          const& Ci = 1,
											    double                                          const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched TMD PDFs in b-space. It
   * assumes the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD PDFs
   */
  std::function<Set<Distribution>(double const&)> MatchTmdPDFs(std::map<int,TmdObjects>                        const& TmdObj,
							       std::function<Set<Distribution>(double const&)> const& CollPDFs,
							       std::function<double(double const&)>            const& Alphas,
							       int                                             const& PerturbativeOrder,
							       double                                          const& Ci = 1);

  /**
   * @brief Function that returns the matched TMD FFs in b-space. It
   * assumes the zeta-prescription.
   * @param TmdObj: the TMD objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD FFs
   */
  std::function<Set<Distribution>(double const&)> MatchTmdFFs(std::map<int,TmdObjects>                        const& TmdObj,
							      std::function<Set<Distribution>(double const&)> const& CollPDFs,
							      std::function<double(double const&)>            const& Alphas,
							      int                                             const& PerturbativeOrder,
							      double                                          const& Ci = 1);

  /**
   * @brief Function that returns the evolution factors for gluon and quarks.
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return std::vector<double>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. The 0-th component contains the
   * gluon evolution vactor, the remaining 12, from 1 to 12, are all
   * equal and represent the quark evolution factors.
   */
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactors(std::map<int,TmdObjects>             const& TmdObj,
												   std::function<double(double const&)> const& Alphas,
												   int                                  const& PerturbativeOrder,
												   double                               const& Ci = 1,
												   double                               const& IntEps = 1e-7);
  ///@}
}
