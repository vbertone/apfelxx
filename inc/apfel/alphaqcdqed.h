//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"
#include "apfel/matrix.h"
#include "apfel/tools.h"

#include <functional>

namespace apfel
{
  /**
   * @brief The AlphaQCDQED is a specialization class of the
   * MatchedEvolution class for the computation of the mixed evolution
   * of QCD and QED.
   */
  class AlphaQCDQED: public MatchedEvolution<matrix<double>>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQCDQED() = delete;

    /**
     * @brief AlphaQCDQED constructor.
     * @param AlphaQCDRef: the reference value of the QCD coupling
     * @param AlphaQEDRef: the reference value of the QED coupling
     * @param MuRef: the reference value of the scale
     * @param QuarkThresholds: vector of quark thresholds
     * @param LepThresholds: vector of lepton thresholds
     * @param pt: perturbative order
     * @param nsteps: number of steps of the ODE solver
     * @note No displaced thresholds allowed.
     */
    AlphaQCDQED(double              const& AlphaQCDRef,
                double              const& AlphaQEDRef,
                double              const& MuRef,
                std::vector<double> const& QuarkThresholds,
                std::vector<double> const& LeptThresholds,
                int                 const& pt,
                int                 const& nsteps = 10);
    ///@}

    /**
     * @brief Function for the computation of the matching.
     * @param Up: tells whether the matching is upward or not (downward)
     * @param nfl: total number of active flavours including both quarks and leptons
     * @param Coup: value of the coupling to be matched
     * @return The matched value of the strong coupling \f$\alpha_s\f$ at the threshold
     * @note This is a dummy function required only to instantiate the
     * pure-virtual 'MatchObject' of the mother class
     * 'MatchedEvolution': it always returns the input coupling. This
     * means that the evolution is assumed to be continuos at the
     * thresholds.
     */
    matrix<double> MatchObject(bool const& Up, int const& nfl, matrix<double> const& Coup) const;

    /**
     * @brief Function that returns QCDxQED \f$\beta\f$ function matrix.
     * @param nfl: total number of active flavours including both quarks and leptons
     * @param as: value of the couplings
     * @return The the value of the QCDxQED \f$\beta\f$ function matrix
     */
    matrix<double> Derivative(int const& nfl, double const&, matrix<double> const& as) const;

    /**
     * @brief Function for the computation of the single coefficients
     * of the expansion of the QCDxQED \f$\beta\f$ function.
     * @param pt: perturbative order
     * @param nf: number of active quark flavours
     * @param nl: number of active lepton flavours
     * @return The pt-th coefficient of the QCDxQED \f$\beta\f$ function.
     */
    matrix<double> betaQCDQED(int const& pt, int const& nf, int const& nl) const;

  private:
    int                                                              const _pt;                    //!< Perturbative order
    std::function<matrix<double>(int const&, matrix<double> const&)>       _BetaFunction;          //!< Beta function
  };
}
