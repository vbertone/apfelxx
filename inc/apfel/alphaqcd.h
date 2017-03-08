//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/matchedevolution.h"

#include <array>

using std::array;

namespace apfel
{
  /**
   * @brief The AlphaQCD class.
   *
   * A specialization class of the MatchedEvolution class
   * for the computation of the QCD coupling running.
   */
  class AlphaQCD: public MatchedEvolution<double>
  {
  public:

    AlphaQCD() =  delete;

    /**
     * @brief AlphaQCD default constructor.
     * @param AlphaRef the reference value of the coupling.
     * @param MuRef the reference value of the scale.
     * @param Masses vector of masses.
     * @param Thresholds vector of thresholds.
     * @param nsteps number of steps of the ODE solver.
     */
    AlphaQCD(double         const& AlphaRef,
	     double         const& MuRef,
	     vector<double> const& Masses,
	     vector<double> const& Thresholds,
	     int            const& pt,
	     int            const& nstep = 10);

    /**
     * @brief AlphaQCD default constructor.
     * @param AlphaRef the reference value of the coupling.
     * @param MuRef the reference value of the scale.
     * @param Masses vector of masses.
     * @param nsteps number of steps of the ODE solver.
     */
    AlphaQCD(double         const& AlphaRef,
	     double         const& MuRef,
	     vector<double> const& Masses,
	     int            const& pt,
	     int            const& nsteps = 10);

    /**
     * @brief Function for the computation of the coupling given nf. This function can be overriden.
     * @param nf number of active flavours.
     * @param as0 starting value of the coupling.
     * @param mu02 initial squared scale.
     * @param mu2 final squared scale.
     * @return value of the coupling at mu2.
     */
    double EvolveObject(int const& nf, double const& mu02, double const& mu2, double const& as0) const;

    /**
     * @brief Function for the computation of the matching. This function can be overriden.
     * @param Up tells whether the matching is upward or not (downward).
     * @param nf number of active flavours.
     * @param Coup value of the coupling to be matched.
     * @return the matched value of the coupling.
     */
    double MatchObject(bool const& Up, int const& nf, double const& Coup) const;

    /**
     * @brief Function for the computation of the terms of the QCD beta function.
     * @param pt perturnative order.
     * @param nf number of active flavours.
     * @return the pt-th coefficient of the beta function.
     */
    double betaQCD(int const& pt, int const& nf) const;

    /**
     * @brief Function for the computation of the full QCD beta function.
     * @param as value of the coupling.
     * @param nf number of active flavours.
     * @return the value of the beta function.
     */
    double fbeta(double const& as, array<double,3> const& bQCD) const;

    /**
     * @brief Function that returns the perturbative order.
     */
    int const& GetPerturbativeOrder() const { return _pt; }

    /**
     * @brief Function that returns the number of steps.
     */
    int const& GetNumberOfSteps()     const { return _nstep; }

  private:
    int                      const _pt;
    int                      const _nstep;
    array<array<double,3>,4>       _bQCD;
  };
}
