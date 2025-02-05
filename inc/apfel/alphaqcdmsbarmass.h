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
   * @brief The AlphaQCDMSbarMass is a specialization class of the
   * MatchedEvolution class for the computation of the QCD coupling
   * running with MSbar (rather than pole) heavy-quark thresholds.
   */
  class AlphaQCDMSbarMass: public MatchedEvolution<double>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQCDMSbarMass() = delete;

    /**
     * @brief AlphaQCDMSbarMass constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param Thresholds: vector of thresholds
     * @param pt: perturbative order
     * @param nsteps: number of steps of the ODE solver
     */
    AlphaQCDMSbarMass(double              const& AlphaRef,
                      double              const& MuRef,
                      std::vector<double> const& Masses,
                      std::vector<double> const& Thresholds,
                      int                 const& pt,
                      int                 const& nsteps = 10);

    /**
     * @brief AlphaQCDMSbarMass constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param pt: perturbative order
     * @param nsteps: number of steps of the ODE solver
     * @note This constructor assumes that masses and thresholds coincide.
     */
    AlphaQCDMSbarMass(double              const& AlphaRef,
                      double              const& MuRef,
                      std::vector<double> const& Masses,
                      int                 const& pt,
                      int                 const& nsteps = 10);
    ///@}

    /**
     * @brief Function for the computation of the matching.
     * @param Up: tells whether the matching is upward or not (downward)
     * @param nf: number of active flavours
     * @param Coup: value of the coupling to be matched
     * @return The matched value of the strong coupling \f$\alpha_s\f$ at the threshold
     */
    double MatchObject(bool const& Up, int const& nf, double const& Coup) const;

    /**
     * @brief Function that returns QCD \f$\beta\f$ function.
     * @param nf: number of active flavours
     * @param as: value of the coupling
     * @return The the value of the QCD \f$\beta\f$ function
     */
    double Derivative(int const& nf, double const&, double const& as) const;

    /**
     * @brief Function for the computation of the single coefficients
     * of the expansion of the QCD \f$\beta\f$ function.
     * @param pt: perturbative order
     * @param nf: number of active flavours
     * @return The pt-th coefficient of the QCD \f$\beta\f$ function.
     */
    double betaQCD(int const& pt, int const& nf) const;

  private:
    int                                                           const _pt;                    //!< Perturbative order
    std::function<double(bool const&, int const&, double const&)>       _MatchingConditions;    //!< Matching condition functions
    std::function<double(int const&, double const&)>                    _BetaFunction;          //!< Beta function
  };
}
