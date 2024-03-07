//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/operator.h"
#include "apfel/set.h"
#include "apfel/dglapbuilder.h"
#include "apfel/tabulateobject.h"
#include "apfel/constants.h"

namespace apfel
{
  struct GtmdObjects
  {
    double                                    Threshold;
    double                                    xi;
    std::map<int, double>                     Beta;
    std::map<int, double>                     GammaFq;
    std::map<int, double>                     GammaFg;
    std::map<int, double>                     GammaK;
    std::map<int, std::vector<double>>        KCS;
    std::map<int, std::vector<Set<Operator>>> MatchingFunctions;
  };

  /**
   * @name GTMD object initializers
   * Collection of functions that initialise GtmdObjects structure for
   * the perturbartive evolution and matching currently available.
   */
  ///@{
  /**
   * @brief The InitializeGtmdObjects function precomputes the
   * perturbative coefficients required for the evolution and matching
   * of GTMD and store them into a 'GtmdObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of GtmdObjects, one for each possible nf
   */
  std::map<int, GtmdObjects> InitializeGtmdObjects(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& xi,
                                                   double              const& IntEps = 1e-5);
  ///@}

  /**
   * @name GTMD builders
   * Collection of functions that build a GTMD distributions as
   * Set<Distribution>-valued functions. These functions perform
   * evolution and matching either separately or alltogether.
   */
  ///@{
  /**
   * @brief Function that returns the matched and evolved GTMDs in
   * b-space as functions of the final scale and rapidity.
   * @param GtmdObj: the GTMD objects
   * @param CollGPDs: the set of collinear GPDs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved GTMDs.
   */
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildGtmds(std::map<int, GtmdObjects>                      const& GtmdObj,
                                                                                           std::function<Set<Distribution>(double const&)> const& CollGPDs,
                                                                                           std::function<double(double const&)>            const& Alphas,
                                                                                           int                                             const& PerturbativeOrder,
                                                                                           double                                          const& Ci = 1,
                                                                                           double                                          const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched TMD GPDs in b-space.
   * @param GtmdObj: the GTMD objects
   * @param CollGPDs: the set of collinear GPDs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched GTMDs.
   */
  std::function<Set<Distribution>(double const&)> MatchGtmds(std::map<int, GtmdObjects>                      const& GtmdObj,
                                                             std::function<Set<Distribution>(double const&)> const& CollGPDs,
                                                             std::function<double(double const&)>            const& Alphas,
                                                             int                                             const& PerturbativeOrder,
                                                             double                                          const& Ci = 1);

  /**
   * @brief Function that returns the mathing functions for the GTMDs.
   * @param GtmdObj: the GTMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Operator>-valued function of the scale mu
   * corresponding to the set of matching functions for GPDs in the
   * evolution basis.
   */
  std::function<Set<Operator>(double const&)> MatchingFunctions(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                std::function<double(double const&)> const& Alphas,
                                                                int                                  const& PerturbativeOrder,
                                                                double                               const& Ci = 1);
  /// @cond UNNECESSARY
  /**
   * @brief Function that returns the evolution factors for gluon and quarks.
   * @param GtmdObj: the GTMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return std::vector<double>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. The 0-th component contains the
   * gluon evolution factor, the remaining 12, from 1 to 12, are all
   * equal and represent the quark evolution factors.
   */
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactors(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                                                   std::function<double(double const&)> const& Alphas,
                                                                                                   int                                  const& PerturbativeOrder,
                                                                                                   double                               const& Ci = 1,
                                                                                                   double                               const& IntEps = 1e-7);

  /**
   * @brief Function that returns the evolution factor for quarks.
   * @param GtmdObj: the GTMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. It returns the quark evolution
   * factor.
   */
  std::function<double(double const&, double const&, double const&)> QuarkEvolutionFactor(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                                          std::function<double(double const&)> const& Alphas,
                                                                                          int                                  const& PerturbativeOrder,
                                                                                          double                               const& Ci = 1,
                                                                                          double                               const& IntEps = 1e-7);

  /**
   * @brief Function that returns the evolution factor for the gluon.
   * @param GtmdObj: the GTMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. It returns the gluon evolution
   * factor.
   */
  std::function<double(double const&, double const&, double const&)> GluonEvolutionFactor(std::map<int, GtmdObjects>           const& GtmdObj,
                                                                                          std::function<double(double const&)> const& Alphas,
                                                                                          int                                  const& PerturbativeOrder,
                                                                                          double                               const& Ci = 1,
                                                                                          double                               const& IntEps = 1e-7);
  /// @endcond
  ///@}
}
