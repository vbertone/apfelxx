//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <cmath>
#include <vector>

using std::vector;

namespace apfel
{

  /**
   * @brief The MatchedEvolution class is a template mother class for
   * the computation of the running of a generic quantity in a
   * VFNS. It provides the basic ingredients for the computation and
   * the heavy-quark threshold matching of the running of a given
   * object.
   */
  template<class T>
  class MatchedEvolution
  {
  public:

    MatchedEvolution() = delete;

    /**
     * @brief The MatchedEvolution default constructor.
     * @param ObjRef: the reference object
     * @param MuRef: the reference scale
     * @param Thresholds: vector with the heavy quark threholds
     * @param nsteps: number of steps of the ODE solver (default: 10)
     */
    MatchedEvolution(T              const& ObjRef,
		     double         const& MuRef,
		     vector<double> const& Thresholds,
		     int            const& nsteps = 10);

    /**
     * @brief Virtual function for the computation of the evolution.
     * @param nf: the number of active flavours
     * @param Obj0: the starting object
     * @param mu02: the squared starting scale
     * @param mu2: the squared final scale
     * @return the object evolved at the scale mu2
     */
    virtual T EvolveObject(int const& nf, double const& mu02, double const& mu2, T const& Obj0) const;

    /**
     * @brief Pure virtual function for the computation of the matching.
     * @param Up: the direction of the matching: "true" = upward, "false" = downward
     * @param nf: the number of active flavours
     * @param Obj: the object to be matched
     * @return the matched object on the other side of the threshold
     */
    virtual T MatchObject(bool const& Up, int const& nf, T const& Obj) const = 0;

    /**
     * @brief Pure virtual function for the computation of the r.h.s. of a ODE (the derivative).
     * @param nf: the number of active flavours
     * @param Mu: the scale at which the derivative is computed
     * @param Obj: the object used to compute the derivative
     * @return the r.h.s. of the ODE
     */
    virtual T Derivative(int const& nf, double const& Mu, T const& Obj) const = 0;

    /**
     * @brief Function that returns the evolved object.
     * @param mu: the final scale
     * @return the evolved object.
     */
    T Evaluate(double const& mu) const;

    /**
     * @name Getters
     */
    ///@{
    /**
     * @brief Function that returns the reference value of the object
     */
    T const& GetObjectRef() const { return _ObjRef; }

    /**
     * @brief Function that returns the reference scale
     */
    double const& GetMuRef() const { return _MuRef; }

    /**
     * @brief Function that returns the values of the thresholds.
     */
    vector<double> const& GetThresholds() const { return _Thresholds; }

    /**
     * @brief Function that returns the number of steps.
     */
    int const& GetNumberOfSteps() const { return _nsteps; }
    ///@}

    /**
     * @name Setters
     */
    ///@{
    /**
     * @brief Function that sets the reference value of the object
     * @param ObjRef: the reference object
     */
    void SetObjectRef(T const& ObjRef) { _ObjRef = ObjRef; }

    /**
     * @brief Function that sets the reference scale
     * @param MuRef: the reference scale
     */
    void SetMuRef(double const& MuRef) { _MuRef2 = MuRef * MuRef; _LogMuRef2 = log(_MuRef2); }

    /**
     * @brief Function that sets the number of steps of the RK algorithm.
     * @param nsteps: the number of steps
     */
    void SetNumberOfSteps(int const& nsteps) { _nsteps = nsteps; }
    ///@}
  protected:
    T              _ObjRef;         //<! Reference value of the object
    double         _MuRef;          //<! Reference scale of the object
    double         _MuRef2;         //<! Squared reference scale of the object
    double         _LogMuRef2;      //<! Log of the squared reference scale of the object
    vector<double> _Thresholds;     //<! Values of the thresholds
    int            _nsteps;         //<! Number of steps of the RK algorithm
    vector<double> _Thresholds2;    //<! Squared quark threholds
    vector<double> _LogThresholds2; //<! Log of the squared quark threholds
  };
}
