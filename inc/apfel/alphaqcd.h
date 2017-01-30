//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <array>
#include <apfel/matchedevolution.h>
using std::array;

namespace apfel
{
  /**
   * @brief The AlphaQCD class.
   *
   * A specialization class of the Coupling class
   * for the computation of the QCD coupling.
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
     */
    AlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, int const& pt, int const& nstep = 10);

    /**
     * @brief Function for the computation of the coupling given nf. This function can be overriden.
     * @param nf number of active flavours.
     * @param as0 starting value of the coupling.
     * @param mu02 initial squared scale.
     * @param mu2 final squared scale.
     * @return value of the coupling at mu2.
     */
    double EvolveObject(int const& nf, double const& as0, double const& mu02, double const& mu2) const;

    /**
     * @brief Function for the computation of the matching. This function can be overriden.
     * @param Up tells whether the matching is upward or not (downward).
     * @param Coup value of the coupling to be matched.
     * @param LogKth value of ln( muth2 / m2 ).
     * @return the matched value of the coupling.
     */
    double MatchObject(bool const& Up, double const& Coup, double const& LogKth) const;

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

  private:
    const int _pt;
    const int _nstep;
  };
}
