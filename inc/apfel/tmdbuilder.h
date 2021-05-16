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

namespace apfel
{
  /**
   * @brief Enumerator for the logarithmic accuracy
   */
  enum LogAccuracy: int {NNNLLp = -3, NNLLp = -2, NLLp = -1, LL = 0, NLL = 1, NNLL = 2, NNNLL = 3};

  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to perform the TMD evolution, matching to the collinear
   * PDFs, and computation of cross sections, i.e. perturbative
   * coefficients of matching functions, all anomalous dimensions, and
   * hard functions.
   */
  struct TmdObjects
  {
    double                                       Threshold;
    std::map<int, double>                        Beta;
    std::map<int, double>                        GammaFq;
    std::map<int, double>                        GammaFg;
    std::map<int, double>                        GammaK;
    std::map<int, std::vector<double>>           KCS;
    std::map<int, std::vector<Set<Operator>>>    MatchingFunctionsPDFs;
    std::map<int, std::vector<Set<Operator>>>    MatchingFunctionsFFs;
    std::map<std::string, std::map<int, double>> HardFactors;
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
   * @param nnnlo: switch to compute the NNNLO matching functions (default: true)
   * @return A map DglapObject objects, one for each possible nf
   */
  std::map<int, TmdObjects> InitializeTmdObjects(Grid                const& g,
                                                 std::vector<double> const& Thresholds,
                                                 double              const& IntEps = 1e-5,
                                                 bool                const& nnnlo = true);

  /**
   * @brief As the InitializeTmdObjects but with the computation of
   * the NNNLO corrections to the matching conditions set to
   * zero. This saves much computation time usually unneded.
   */
  std::map<int, TmdObjects> InitializeTmdObjectsLite(Grid                const& g,
                                                     std::vector<double> const& Thresholds,
                                                     double              const& IntEps = 1e-5);
  ///@}

  /**
   * @name TMD builders
   * Collection of functions that build a TMD distributions (both PDFs
   * and FF) as Set<Distribution>-valued functions. These functions
   * perform evolution and matching either separately or alltogether.
   * Also a function for the computation of the hard factors is
   * provided.
   */
  ///@{
  /**
   * @brief Function that returns the matched and evolved TMD PDFs in
   * b-space as functions of the final scale and rapidity.
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
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdPDFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                                                             std::function<Set<Distribution>(double const&)> const& CollPDFs,
                                                                                             std::function<double(double const&)>            const& Alphas,
                                                                                             int                                             const& PerturbativeOrder,
                                                                                             double                                          const& Ci = 1,
                                                                                             double                                          const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched and evolved TMD FFs in
   * b-space as functions of the final scale and rapidity.
   * @param TmdObj: the TMD objects
   * @param CollFFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved TMD FFs
   */
  std::function<Set<Distribution>(double const&, double const&, double const&)> BuildTmdFFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                                                            std::function<Set<Distribution>(double const&)> const& CollFFs,
                                                                                            std::function<double(double const&)>            const& Alphas,
                                                                                            int                                             const& PerturbativeOrder,
                                                                                            double                                          const& Ci = 1,
                                                                                            double                                          const& IntEps = 1e-7);

  /**
   * @brief Function that returns the matched TMD PDFs in b-space.
   * @param TmdObj: the TMD objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD PDFs
   */
  std::function<Set<Distribution>(double const&)> MatchTmdPDFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                               std::function<Set<Distribution>(double const&)> const& CollPDFs,
                                                               std::function<double(double const&)>            const& Alphas,
                                                               int                                             const& PerturbativeOrder,
                                                               double                                          const& Ci = 1);

  /**
   * @brief Function that returns the matched TMD FFs in b-space.
   * @param TmdObj: the TMD objects
   * @param CollPDFs: the set of collinear PDFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD FFs
   */
  std::function<Set<Distribution>(double const&)> MatchTmdFFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                              std::function<Set<Distribution>(double const&)> const& CollPDFs,
                                                              std::function<double(double const&)>            const& Alphas,
                                                              int                                             const& PerturbativeOrder,
                                                              double                                          const& Ci = 1);

  /**
   * @brief Function that returns the mathing functions for the TMD PDFs.
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Operator>-valued function of the scale mu
   * corresponding to the set of matching functions for PDFs in the
   * evolution basis.
   */
  std::function<Set<Operator>(double const&)> MatchingFunctionsPDFs(std::map<int, TmdObjects>            const& TmdObj,
                                                                    std::function<double(double const&)> const& Alphas,
                                                                    int                                  const& PerturbativeOrder,
                                                                    double                               const& Ci = 1);

  /**
   * @brief Function that returns the mathing functions for the TMD FFs.
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Operator>-valued function of the scale mu
   * corresponding to the set of matching functions for FFs in the
   * evolution basis.
   */
  std::function<Set<Operator>(double const&)> MatchingFunctionsFFs(std::map<int, TmdObjects>            const& TmdObj,
                                                                   std::function<double(double const&)> const& Alphas,
                                                                   int                                  const& PerturbativeOrder,
                                                                   double                               const& Ci = 1);

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
   * gluon evolution factor, the remaining 12, from 1 to 12, are all
   * equal and represent the quark evolution factors.
   */
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactors(std::map<int, TmdObjects>            const& TmdObj,
                                                                                                   std::function<double(double const&)> const& Alphas,
                                                                                                   int                                  const& PerturbativeOrder,
                                                                                                   double                               const& Ci = 1,
                                                                                                   double                               const& IntEps = 1e-7);

  /**
   * @brief Function that returns the evolution factor for quarks.
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. It returns the quark evolution
   * factor.
   */
  std::function<double(double const&, double const&, double const&)> QuarkEvolutionFactor(std::map<int, TmdObjects>            const& TmdObj,
                                                                                          std::function<double(double const&)> const& Alphas,
                                                                                          int                                  const& PerturbativeOrder,
                                                                                          double                               const& Ci = 1,
                                                                                          double                               const& IntEps = 1e-7);

  /**
   * @brief Function that returns the evolution factor for the gluon.
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta;. It returns the gluon evolution
   * factor.
   */
  std::function<double(double const&, double const&, double const&)> GluonEvolutionFactor(std::map<int, TmdObjects>            const& TmdObj,
                                                                                          std::function<double(double const&)> const& Alphas,
                                                                                          int                                  const& PerturbativeOrder,
                                                                                          double                               const& Ci = 1,
                                                                                          double                               const& IntEps = 1e-7);

  /**
   * @brief Function that returns the hard factor.
   * @param Process: the string corresponding to the process requested
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Cf: the final scale-variation factor (default: 1)
   * @return double-valued function of the final renormalisation scale &mu;
   */
  std::function<double(double const&)> HardFactor(std::string                          const& Process,
                                                  std::map<int, TmdObjects>            const& TmdObj,
                                                  std::function<double(double const&)> const& Alphas,
                                                  int                                  const& PerturbativeOrder,
                                                  double                               const& Cf = 1);
  ///@}

  /**
   * @name Hard factors
   * Collection of hard functions used in the defintion of the
   * qT-dependent cross sections in TMD factorisation.
   */
  ///@{
  /**
   * @brief Perturbative hard factor for Drell-Yan.
   * @param PerturbativeOrder: the perturbative order
   * @param Alphas: the value of the strong coupling
   * @param nf: the number of active flavours
   * @param kappa: the ration between hard scale and renormalusation scale
   * @return The hard factor for Drell-Yan.
   */
  double HardFactorDY(int const& PerturbativeOrder, double const& Alphas, int const& nf, double const& kappa);
  ///@}
}
