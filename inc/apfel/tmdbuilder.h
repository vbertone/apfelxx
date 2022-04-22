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
   * @brief Enumerator for the jet algoritms fot the jet TMDs
   */
  enum JetAlgorithm: int {CONE = 0, KT = 1};

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
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @param nnnlo: switch to compute the NNNLO matching functions (default: false)
   * @return A map of TmdObject objects, one for each possible nf
   */
  std::map<int, TmdObjects> InitializeTmdObjects(Grid                const& g,
                                                 std::vector<double> const& Thresholds,
                                                 double              const& IntEps = 1e-5,
                                                 bool                const& nnnlo = false);

  /**
   * @brief The InitializeTmdObjectsDYResScheme function precomputes
   * the perturbative coefficients required for the evolution and
   * matching of TMD PDFs and FFs and store them into a 'TmdObjects'
   * structure. This function applies a resummation-scheme
   * transformation to produce the scheme often used in qT resummation
   * that has H = 1.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of TmdObject objects, one for each possible nf
   */
  std::map<int, TmdObjects> InitializeTmdObjectsDYResScheme(Grid                const& g,
                                                            std::vector<double> const& Thresholds,
                                                            double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeTmdObjectsBM function precomputes the
   * perturbative coefficients required for the evolution and matching
   * of the (gluon) Boer-Mulders TMD PDF and store them into a
   * 'TmdObjects' structure. For now, quark and FF TMDs are not filled
   * in.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return A map of TmdObject objects, one for each possible nf
   */
  std::map<int, TmdObjects> InitializeTmdObjectsBM(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeTmdObjectsSivers function precomputes the
   * perturbative coefficients required for the evolution and matching
   * of the quark Sivers TMD PDF and store them into a 'TmdObjects'
   * structure. For now, gluon and FF TMDs (i.e. the Collins TMDs) are
   * not filled in. In addition, the matching is only present up to
   * one loop.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return A map of TmdObject objects, one for each possible nf
   */
  std::map<int, TmdObjects> InitializeTmdObjectsSivers(Grid                const& g,
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
   * @param Ci: the initial-scale variation factor (default: 1)
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
   * @param Ci: the initial-scale variation factor (default: 1)
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
   * @brief Function that returns the TMD of a jet in b-space as
   * functions of the final scale and rapidity.
   * @param TmdObj: the TMD objects
   * @param JetAlgo: the jet algorithm to be used
   * @param JetR: the jet radius
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param CJ: jet-scale variation factor (default: 1)
   * @param Ci: the initial-scale variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB>, the final renormalisation scale &mu;, and the
   * final rapidity scale &zeta; representing the evolved jet TMD
   */
  std::function<double(double const&, double const&, double const&)> BuildTmdJet(std::map<int, TmdObjects>            const& TmdObj,
                                                                                 JetAlgorithm                         const& JetAlgo,
                                                                                 double                               const& JetR,
                                                                                 std::function<double(double const&)> const& Alphas,
                                                                                 int                                  const& PerturbativeOrder,
                                                                                 double                               const& CJ = 1,
                                                                                 double                               const& Ci = 1,
                                                                                 double                               const& IntEps = 1e-7);

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
   * @param CollFFs: the set of collinear FFs to be matched
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial-scale variation factor
   * @return Set<Distribution>-valued function of the impact parameter
   * b<SUB>T</SUB> representing the matched TMD FFs
   */
  std::function<Set<Distribution>(double const&)> MatchTmdFFs(std::map<int, TmdObjects>                       const& TmdObj,
                                                              std::function<Set<Distribution>(double const&)> const& CollFFs,
                                                              std::function<double(double const&)>            const& Alphas,
                                                              int                                             const& PerturbativeOrder,
                                                              double                                          const& Ci = 1);

  /**
   * @brief Function that returns the jet TMD in b-space at the
   * initial scale.
   * @param TmdObj: the TMD objects
   * @param JetAlgo: the jet algorithm
   * @param tR: tangent of half the jet radius (tan(R/2))
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param CJ: jet-scale variation factor (default: 1)
   * @param Ci: the initial-scale variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB> representing the low-scale jet TMD
   */
  std::function<double(double const&, double const&)> MatchTmdJet(std::map<int, TmdObjects>            const& TmdObj,
                                                                  JetAlgorithm                         const& JetAlgo,
                                                                  double                               const& tR,
                                                                  std::function<double(double const&)> const& Alphas,
                                                                  int                                  const& PerturbativeOrder,
                                                                  double                               const& CJ = 1,
                                                                  double                               const& Ci = 1,
                                                                  double                               const& IntEps = 1e-7);

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
   * @brief Function that returns the evolution factors for gluon and
   * quarks. As compared to "EvolutionFactors", this function isolates
   * the double logs into gammaK. This is reminiscent of the
   * qT-resummation typical way of computing the Sudakov form factor.
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
  std::function<std::vector<double>(double const&, double const&, double const&)> EvolutionFactorsK(std::map<int, TmdObjects>            const& TmdObj,
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
   * @brief Function that returns the perturbative part of the
   * Collins-Soper kernel for quarks.
   * @param TmdObj: the TMD objects
   * @param Alphas: the strong coupling function
   * @param PerturbativeOrder: the perturbative order
   * @param Ci: the initial scale-variation factor (default: 1)
   * @param IntEps: the integration accuracy (default: 10<SUP>-7</SUP>)
   * @return double-valued function of the impact parameter
   * b<SUB>T</SUB> and of the the final renormalisation scale &mu;. It
   * returns perturbative part of the Collis-Soper kernel for quarks.
   */
  std::function<double(double const&, double const&)> CollinsSoperKernel(std::map<int, TmdObjects>            const& TmdObj,
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
   * @note TO BE REMOVED!
   */
  double HardFactorDY(int const& PerturbativeOrder, double const& Alphas, int const& nf, double const& kappa);
  ///@}
}
