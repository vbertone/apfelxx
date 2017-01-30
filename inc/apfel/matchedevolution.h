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
   * @brief The MatchedEvolution class for the computation of the running of a generic quantity. It can be couplings as well as distributions.
   *
   * This mother class provides the basic ingredients for the computation
   * and the threshold matching of the running of a given object.
   */
  template<class T>
  class MatchedEvolution
  {
  public:

    MatchedEvolution() = delete;

    /**
     * @brief The default constructor that takes the reference object and the reference scale.
     * @param ObjRef reference object
     * @param MuRef reference scale
     * @param Masses vector with the heavy quark masses
     * @param Thresholds vector with the heavy quark threholds
     */
    MatchedEvolution(T const& ObjRef, double const& MuRef, vector<double> const& Masses, vector<double> const& Thresholds);

    /**
     * @brief The default constructor that takes the reference value of the object and the reference scale (assumes equal masses and thresholds).
     * @param ObjRef reference object
     * @param MuRef reference scale
     * @param Masses vector with the heavy quark masses
     */
    MatchedEvolution(T const& ObjRef, double const& MuRef, vector<double> const& Masses);

    /**
     * @brief Virtual function for the computation of the evolution of the object with nf flavours.
     * @param nf number of flavours (constant during the evolution).
     * @param Obj0 starting object.
     * @param mu02 squared starting scale.
     * @param mu2 squared final scale.
     * @return the object at the scale mu2.
     */
    virtual T EvolveObject(int const& nf, T const& Obj0, double const& mu02, double const& mu2) const = 0;

    /**
     * @brief Virtual function for the computation of the matching.
     * @param Up direction of the matching "true" = upward, "false" = downward
     * @param Obj object to be matched
     * @param LogKth value of ln(muth2/m2), where muth2 is the threshold and m2 the mass, both squared  
     * @return the matched object.
     */
    virtual T MatchObject(bool const& Up, T const& Obj, double const& LogKth) const = 0;

    /**
     * @brief Function that returns the evolved Object.
     * @param mu final scale
     * @return the evolved Object.
     */
    T GetObject(double const& mu) const;

    /**
     * @brief Function that returns the values of the thresholds.
     */
    vector<double> const& GetThresholds() const { return _Thresholds; }

    /**
     * @brief Function that returns the values of the masses.
     */
    vector<double> const& GetMasses() const { return _Masses; }

    /**
     * @brief Function that sets the reference value of the object
     * @param ObjRef
     */
    void SetObjectRef(T const& ObjRef) { _ObjRef = ObjRef; }

    /**
     * @brief Function that sets the reference scale
     * @param MuRef
     */
    void SetMuRef(double const& MuRef) { _MuRef2 = MuRef * MuRef; }

  protected:
    T              _ObjRef;       //<! Reference value of the object
    double         _MuRef2;       //<! Squared reference scale of the object
    vector<double> _Masses;       //<! Values of the masses
    vector<double> _Thresholds;   //<! Values of the thresholds
    vector<double> _Thresholds2;  //<! Squared quark threholds
    vector<double> _LogTh2M2;     //<! Log of the squared threholds over squared masses

  };
}
