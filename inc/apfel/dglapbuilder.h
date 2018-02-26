//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/evolutionbasisqcd.h"
#include "apfel/matchingbasisqcd.h"

#include <functional>
#include <vector>

using std::function;
using std::vector;

namespace apfel
{
  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to perform the DGLAP evolution, i.e. perturbative
   * coefficients of splitting functions and matching conditions, and
   * the heavy quark thresholds.
   */
  struct DglapObjects
  {
    double Threshold;
    map<int,Set<Operator>> SplittingFunctions;
    map<int,Set<Operator>> MatchingConditions;
  };

  /**
   * @name DGLAP object initializers
   * Collection of functions that initialise DglapObjects structure
   * for the different kinds of evolution currently available.
   */
  ///@{
  /**
   * @brief The InitializeDglapObjectsQCD function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   */
  map<int,DglapObjects> InitializeDglapObjectsQCD(Grid           const& g,
						  vector<double> const& Masses,
						  vector<double> const& Thresholds,
						  bool           const& OpEvol = false,
						  double         const& IntEps = 1e-5);
  /**
   * @brief The InitializeDglapObjectsQCD function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  map<int,DglapObjects> InitializeDglapObjectsQCD(Grid           const& g,
						  vector<double> const& Thresholds,
						  bool           const& OpEvol = false,
						  double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDT function precomputes the
   * perturbative coefficients of time-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   */
  map<int,DglapObjects> InitializeDglapObjectsQCDT(Grid           const& g,
						   vector<double> const& Masses,
						   vector<double> const& Thresholds,
						   bool           const& OpEvol = false,
						   double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDT function precomputes the
   * perturbative coefficients of time-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  map<int,DglapObjects> InitializeDglapObjectsQCDT(Grid           const& g,
						   vector<double> const& Thresholds,
						   bool           const& OpEvol = false,
						   double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of space-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   */
  map<int,DglapObjects> InitializeDglapObjectsQCDtrans(Grid           const& g,
						       vector<double> const& Masses,
						       vector<double> const& Thresholds,
						       bool           const& OpEvol = false,
						       double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of space-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  map<int,DglapObjects> InitializeDglapObjectsQCDtrans(Grid           const& g,
						       vector<double> const& Thresholds,
						       bool           const& OpEvol = false,
						       double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of timelike-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   */
  map<int,DglapObjects> InitializeDglapObjectsQCDTtrans(Grid           const& g,
							vector<double> const& Masses,
							vector<double> const& Thresholds,
							bool           const& OpEvol = false,
							double         const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of time-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param OpEvol: the switsch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10^{-5})
   * @return A map DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  map<int,DglapObjects> InitializeDglapObjectsQCDTtrans(Grid           const& g,
							vector<double> const& Thresholds,
							bool           const& OpEvol = false,
							double         const& IntEps = 1e-5);
  ///@}

  /**
   * @name DGLAP builders
   * Collection of functions that build a Dglap object used to perform
   * the DGLAP evolution of distributions or operators.
   */
  ///@{
  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs the DGLAP evolution for distributions.
   * @param DglapObj: structure with the coefficients of the perturbative objects
   * @param InDistFunc: the distributions at the reference scale
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  unique_ptr<Dglap<Distribution>> BuildDglap(map<int,DglapObjects>                                   const& DglapObj,
					     function<map<int,double>(double const&, double const&)> const& InDistFunc,
					     double                                                  const& MuRef,
					     int                                                     const& PerturbativeOrder,
					     function<double(double const&)>                         const& Alphas,
					     int                                                     const& nsteps = 10);

  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs the DGLAP evolution for operators.
   * @param DglapObj: structure with the coefficients of the perturbative objects
   * @param InDistFunc: the distributions at the reference scale
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  unique_ptr<Dglap<Operator>> BuildDglap(map<int,DglapObjects>           const& DglapObj,
					 double                          const& MuRef,
					 int                             const& PerturbativeOrder,
					 function<double(double const&)> const& Alphas,
					 int                             const& nsteps = 10);

  /**
   * @brief Same as above but the DglapObjects object is a function of
   * a variable. This is useful for applications in which the
   * splitting functions depend on same energy scale.
   */
  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs the DGLAP evolution for distributions.
   * @param DglapObj: DglapObjects-values function that returns the structure with the coefficients of the perturbative objects as function of a scale
   * @param InDistFunc: the distributions at the reference scale
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  unique_ptr<Dglap<Distribution>> BuildDglap(function<DglapObjects(double const&)>                   const& DglapObj,
					     vector<double>                                          const& Thresholds,
					     function<map<int,double>(double const&, double const&)> const& InDistFunc,
					     double                                                  const& MuRef,
					     int                                                     const& PerturbativeOrder,
					     function<double(double const&)>                         const& Alphas,
					     int                                                     const& nsteps = 10);
  ///@}
}
