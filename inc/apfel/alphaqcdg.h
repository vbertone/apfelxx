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
   * @brief The AlphaQCDg is a specialization class of the
   * MatchedEvolution class for the computation of the QCD coupling
   * running using the analytic g functions.
   */
  class AlphaQCDg: public MatchedEvolution<double>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQCDg() = delete;

    /**
     * @brief AlphaQCDg constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param Thresholds: vector of thresholds
     * @param pt: perturbative order
     * @param kappa: resummation scale parameter (default: 1)
     */
    AlphaQCDg(double              const& AlphaRef,
              double              const& MuRef,
              std::vector<double> const& Masses,
              std::vector<double> const& Thresholds,
              int                 const& pt,
              double              const& kappa = 1);

    /**
     * @brief AlphaQCDg constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param pt: perturbative order
     * @param kappa: resummation scale parameter (default: 1)
     * @note This constructor assumes that masses and thresholds coincide.
     */
    AlphaQCDg(double              const& AlphaRef,
              double              const& MuRef,
              std::vector<double> const& Masses,
              int                 const& pt,
              double              const& kappa = 1);
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
     * @brief Dummy function used to overload a purely virtual
     * function. It should be never called.
     */
    double Derivative(int const&, double const&, double const&) const { return 0; };

    /**
     * @brief Function for the computation of the evolution.
     * @param nf: the number of active flavours
     * @param lnmu02: the log of the initial value of the scale
     * @param lnmu2 the log of the final value of the scale
     * @param as0: the starting object
     * @return the object evolved at the scale mu2
     */
    double EvolveObject(int const& nf, double const& lnmu02, double const& lnmu2, double const& as0) const;

  private:
    int                                                           const _pt;                    //!< Perturbative order
    double                                                        const _kappa;                 //!< Resummation-scale parameter
    std::function<double(bool const&, int const&, double const&)>       _MatchingConditions;    //!< Matching condition functions
  };
}
