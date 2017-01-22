//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>
using std::vector;

namespace apfel
{

  /**
   * @brief The Expression class for the computation of the running of the couplings.
   *
   * This mother class provides the basic ingredients for the computation
   * and the threshold matching of the running of a given coupling.
   */
  class Coupling
  {
  public:

    Coupling() = delete;

    /**
     * @brief The default constructor that takes the reference value of the coupling and the reference scale.
     * @param AlphaRef reference value
     * @param MuRef reference scale
     * @param Masses vector with the heavy quark masses
     * @param Thresholds vector with the heavy quark threholds
     */
    Coupling(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, vector<double> const& Thresholds);

    /**
     * @brief The default constructor that takes the reference value of the coupling and the reference scale (assumes equal masses and thresholds).
     * @param AlphaRef reference value
     * @param MuRef reference scale
     * @param Masses vector with the heavy quark masses
     */
    Coupling(double const& AlphaRef, double const& MuRef, vector<double> const& Masses);

    /**
     * @brief Virtual function for the computation of the coupling.
     * @param nf number of flavours (constant during the evolution).
     * @param coupling0 starting value of coupling.
     * @param mu02 squared starting scale.
     * @param mu2 squared final scale.
     * @return the coupling at the scale mu2.
     */
    virtual double Coup(int const& nf, double const& coupling0, double const& mu02, double const& mu2) const = 0;

    /**
     * @brief Virtual function for the computation of the matching.
     * @param Up direction of the matching "true" = upward, "false" = downward
     * @param Coup coupling to be matched
     * @param LogKth value of ln(muth2/m2), where muth2 is the threshold and m2 the mass, both squared  
     * @return the matched coupling.
     */
    virtual double MatchCoupling(bool const& Up, double const& Coup, double const& LogKth) const = 0;

    /**
     * @brief Function that returns the actual coupling.
     * @param mu final scale
     * @return the evolved coupling.
     */
    double GetCoupling(double const& mu) const;

    /**
     * @brief Function that returns the values of the thresholds.
     */
    vector<double> const& GetThresholds() const { return _Thresholds; }

    /**
     * @brief Function that returns the values of the masses.
     */
    vector<double> const& GetMasses() const { return _Masses; }

    /**
     * @brief Function that sets the reference value of the coupling
     * @param AlphaRef
     */
    void SetAlphaRef(double const& AlphaRef) { _AlphaRef = AlphaRef; }

    /**
     * @brief Function that sets the reference scale
     * @param MuRef
     */
    void SetMuRef(double const& MuRef) { _MuRef2 = MuRef * MuRef; }

  protected:
    double         _AlphaRef;     //<! Reference value of the coupling
    double         _MuRef2;       //<! Squared reference scale of the coupling
    vector<double> _Masses;       //<! Values of the masses
    vector<double> _Thresholds;   //<! Values of the thresholds
    vector<double> _Thresholds2;  //<! Squared quark threholds
    vector<double> _LogTh2M2;     //<! Log of the squared threholds over squared masses

  };
}
