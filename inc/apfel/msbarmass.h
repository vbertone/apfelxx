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
   * @brief The MSbarMass is a specialization class of the
   * MatchedEvolution class for the evolution of masses renormalised
   * in the MSbar scheme.
   */
  class MSbarMass: public MatchedEvolution<double>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    MSbarMass() = delete;

    /**
     * @brief MSbarMass constructor.
     * @param MRef: the reference value of the mass
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param pt: perturbative order
     * @param Alphas: the function returning the strong coupling
     * @param nsteps: number of steps of the ODE solver
     * @note Matching conditions assume that the matching thresholds
     * coincide with the RGI values mh(mh).
     */
    MSbarMass(double                               const& MRef,
              double                               const& MuRef,
              std::vector<double>                  const& Masses,
              int                                  const& pt,
              std::function<double(double const&)> const& Alphas,
              int                                  const& nsteps = 10);
    ///@}

    /**
     * @brief Function for the computation of the matching.
     * @param Up: tells whether the matching is upward or not (downward)
     * @param nf: number of active flavours
     * @param Coup: value of the coupling to be matched
     * @return The matched value of the mass at the threshold
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
     * @brief Function for the computation of the \f$\gamma_m\f$
     * anomalous dimension.
     * @param pt: perturbative order
     * @param nf: number of active flavours
     * @return The pt-th coefficient of the QCD \f$\beta\f$ function.
     */
    double AnomalouDimension(int const& pt, int const& nf) const;

  private:
    int                                                           const _pt;                    //!< Perturbative order
    std::function<double(double const&)>                          const _Alphas;                //!< Running coupling
    std::function<double(bool const&, int const&, double const&)>       _MatchingConditions;    //!< Matching condition functions
    std::function<double(int const&, double const&)>                    _AnomalousDimension;    //!< Anomalous dimension
  };
}
