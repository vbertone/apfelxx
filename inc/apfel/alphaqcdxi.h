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
   * @brief The AlphaQCDxi is a specialization class of the
   * MatchedEvolution class for the computation of the QCD coupling
   * running with the possibility ro vary the resummation scale
   * through the parameter xi.
   */
  class AlphaQCDxi: public MatchedEvolution<double>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQCDxi() = delete;

    /**
     * @brief AlphaQCDxi constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param Thresholds: vector of thresholds
     * @param pt: perturbative order
     * @param xi: resummation-scale parameter (default: 1)
     * @param nsteps: number of steps of the ODE solver (default: 10)
     */
    AlphaQCDxi(double              const& AlphaRef,
               double              const& MuRef,
               std::vector<double> const& Masses,
               std::vector<double> const& Thresholds,
               int                 const& pt,
               double              const& xi = 1,
               int                 const& nsteps = 10);

    /**
     * @brief AlphaQCDxi constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param pt: perturbative order
     * @param xi: resummation-scale parameter (default: 1)
     * @param nsteps: number of steps of the ODE solver (default: 10)
     * @note This constructor assumes that masses and thresholds coincide.
     */
    AlphaQCDxi(double              const& AlphaRef,
               double              const& MuRef,
               std::vector<double> const& Masses,
               int                 const& pt,
               double              const& xi = 1,
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
     * @param as: value of the coupling
     * @param nf: number of active flavours
     * @return The the value of the QCD \f$\beta\f$ function
     */
    double Derivative(int const& nf, double const&, double const& as) const;

    /**
     * @brief Function for the computation of the single coefficients of the expansion of the QCD \f$\beta\f$ function.
     * @param pt: perturbative order
     * @param nf: number of active flavours
     * @param lambda: argument of the g-functions
     * @return The pt-th coefficient of the QCD \f$\beta\f$ function.
     */
    double betaQCD(int const& pt, int const& nf, double const& lambda) const;

  private:
    int                                                           const _pt;                    //!< Perturbative order
    double                                                        const _xi;                    //!< Resummation-scale paremeter
    std::function<double(bool const&, int const&, double const&)>       _MatchingConditions;    //!< Matching condition functions
    std::function<double(int const&, double const&)>                    _BetaFunction;          //!< Beta function
  };
}
