//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"
#include "apfel/constants.h"

#include <functional>
#include <complex>

namespace apfel
{
  /**
   * @brief The AlphaQCDLambda is a specialization class of the
   * MatchedEvolution class for the computation of the QCD coupling
   * running in terms of LambdaQCD.
   */
  class AlphaQCDLambda
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQCDLambda() = delete;

    /**
     * @brief AlphaQCDLambda constructor.
     * @param LambdaQCD: the value of LambdaQCD
     * @param nfRef: the number of flavours of LambdaQCD
     * @param Thresholds: vector of thresholds
     * @param pt: perturbative order
     * @param kappa: resummation scale parameter (default: 1)
     */
    AlphaQCDLambda(double              const& LambdaQCD,
                   int                 const& nfRef,
                   std::vector<double> const& Thresholds,
                   int                 const& pt);
    ///@}

    /**
     * @brief Function that returns the evolved object over the full complex plane
     * @param t: the complex final scale expressed as log(m2)
     * @return the evolved object.
     */
    std::complex<double> Evaluate(std::complex<double> const& t) const;

    /**
     * @brief Function that returns the evolved object.
     * @param mu: the final scale
     * @return the evolved object.
     */
    double Evaluate(double const& mu) const;

    /**
     * @brief Function that returns the evolved object using the
     * analytic-perturbation-theory prescription.
     * @param mu: the final scale
     * @param p: the power of the coupling (default: 1)
     * @param tmin: lower integration limit in t (default: -100)
     * @param tmax: upper integration limit in t (default: 30000)
     * @param eps: integration accuracy (default: 10<SUP>-5</SUP>)
     * @return the evolved object.
     */
    double EvaluateAPT(double const& mu, double const& p = 1, double const& tmin = -100, double const& tmax = 30000, double const& eps = eps5) const;

    /**
     * @brief Function that returns the values of LambdaQCD.
     */
    std::vector<double> GetLambdaQCD() const { return _LambdaQCD; };

  private:
    std::vector<double> const _Thresholds; //!< Vector of heavy-quark thresholds
    int                 const _pt;         //!< Perturbative order
    std::vector<double>       _LambdaQCD;  //!< Values of LambdaQCD according to the number of active flavours
    std::vector<double>       _tThresholds; //!< Vector of heavy-quark thresholds in the variable t = ln(mu2)
  };
}
