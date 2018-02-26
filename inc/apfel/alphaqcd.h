//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/matchedevolution.h>
#include <apfel/matrix.h>

#include <functional>

using std::function;

namespace apfel
{
  /**
   * @brief The AlphaQCD class.
   *
   * A specialization class of the MatchedEvolution class for the
   * computation of the QCD coupling running.
   */
  class AlphaQCD: public MatchedEvolution<double>
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    AlphaQCD() = delete;

    /**
     * @brief AlphaQCD default constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param Thresholds: vector of thresholds
     * @param nsteps: number of steps of the ODE solver
     */
    AlphaQCD(double         const& AlphaRef,
	     double         const& MuRef,
	     vector<double> const& Masses,
	     vector<double> const& Thresholds,
	     int            const& pt,
	     int            const& nstep = 10);

    /**
     * @brief AlphaQCD default constructor.
     * @param AlphaRef: the reference value of the coupling
     * @param MuRef: the reference value of the scale
     * @param Masses: vector of masses
     * @param nsteps: number of steps of the ODE solver
     * @note This constructor assumes that masses and thresholds coincide.
     */
    AlphaQCD(double         const& AlphaRef,
	     double         const& MuRef,
	     vector<double> const& Masses,
	     int            const& pt,
	     int            const& nsteps = 10);
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
     * @param as value of the coupling
     * @param nf number of active flavours
     * @return The the value of the QCD \f$\beta\f$ function
     */
    double Derivative(int const& nf, double const&, double const& as) const;

    /**
     * @brief Function for the computation of the single coefficients of the expansion of the QCD \f$\beta\f$ function.
     * @param pt: perturbative order
     * @param nf: number of active flavours
     * @return The pt-th coefficient of the QCD \f$\beta\f$ function.
     */
    double betaQCD(int const& pt, int const& nf) const;

  private:
    int                                                      const _pt;
    matrix<double>                                                 _bQCD;
    function<double(bool const&, int const&, double const&)>       _MatchingConditions;
    function<double(int const&, double const&)>                    _BetaFunction;
  };
}
