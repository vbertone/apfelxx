//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"

#include <functional>

namespace apfel
{
  /**
   * @brief The AlphaQED is a specialization class of the
   * MatchedEvolution class for the computation of the QED coupling
   * running.
   */
  class AlphaQED: public MatchedEvolution<double>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQED() = delete;

    /**
     * @brief AlphaQED constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param QuarkThresholds: vector of quark thresholds
     * @param LeptThresholds: vector of charged-lepton thresholds
     * @param nsteps: number of steps of the ODE solver
     */
    AlphaQED(double              const& AlphaRef,
	     double              const& MuRef,
	     std::vector<double> const& LeptThresholds,
	     std::vector<double> const& QuarkThresholds,
	     int                 const& pt,
	     int                 const& nstep = 10);
    ///@}

    /**
     * @brief Function for the computation of the matching.
     * @param Up: tells whether the matching is upward or not (downward)
     * @param nf: number of active flavours
     * @param Coup: value of the coupling to be matched
     * @return The matched value of the strong coupling \f$\alpha_s\f$ at the threshold
     * @note This is a dummy function required only to instantiate the
     * pure-virtual 'MatchObject' of the mother class
     * 'MatchedEvolution': it always returns the input coupling. This
     * means that the evolution is assumed to be continuos at the
     * thresholds.
     */
    double MatchObject(bool const& Up, int const& nf, double const& Coup) const;

    /**
     * @brief Function that returns QED \f$\beta\f$ function.
     * @param a: value of the coupling
     * @param nfl: total number of active flavours including both quarks and leptons
     * @return The the value of the QED \f$\beta\f$ function
     */
    double Derivative(int const& nfl, double const&, double const& a) const;

    /**
     * @brief Function for the computation of the single coefficients of the expansion of the QED \f$\beta\f$ function.
     * @param pt: perturbative order
     * @param nf: number of active quark flavours
     * @param nl: number of active charged leptons
     * @return The pt-th coefficient of the QED \f$\beta\f$ function.
     */
    double betaQED(int const& pt, int const& nf, int const& nl) const;

  private:
    int                                        const _pt;                    //!< Perturbative order
    std::function<double(int const&, double const&)> _BetaFunction;          //!< Beta function
  };
}
