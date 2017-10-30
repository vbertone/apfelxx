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
   * coefficients of splitting functions and matching conditions.
   */
  struct DglapObjects
  {
    double Threshold;
    map<int,Set<Operator>> SplittingFunctions;
    map<int,Set<Operator>> MatchingConditions;
  };

  /**
   * @brief The InitializeDglapObjectsQCD, precompute the perturbative
   * coefficients of splitting functions and matching conditions and
   * store them in the 'DglapObjects' structure.
   * @param g the grid
   * @param IntEps the integration accuracy
   * @param Masses
   * @param Thresholds
   * @return
   */
  map<int,DglapObjects> InitializeDglapObjectsQCD(Grid           const& g,
						  vector<double> const& Masses,
						  vector<double> const& Thresholds,
						  bool           const& OpEvol = false,
						  double         const& IntEps = 1e-5);

  /**
   * @brief Same as above but assuming that Masses = Thresholds.
   */
  map<int,DglapObjects> InitializeDglapObjectsQCD(Grid           const& g,
						  vector<double> const& Thresholds,
						  bool           const& OpEvol = false,
						  double         const& IntEps = 1e-5);

  // ============================================================================
  /**
   * @brief The BuildDglap, builds the dglap object
   *
   * @param DglapObj structure with the coeffs. of perturbative objects
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param nsteps the number of steps for RK.
   * @return
   */
  unique_ptr<Dglap<Distribution>> BuildDglap(map<int,DglapObjects>                                   const& DglapObj,
					     function<map<int,double>(double const&, double const&)> const& InDistFunc,
					     double                                                  const& MuRef,
					     int                                                     const& PerturbativeOrder,
					     function<double(double const&)>                         const& Alphas,
					     int                                                     const& nsteps = 10);
  /**
   * @brief Same as above but it computes the evolution of the
   * evolution operators.
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
  unique_ptr<Dglap<Distribution>> BuildDglap(function<DglapObjects(double const&)>                   const& DglapObj,
					     vector<double>                                          const& Thresholds,
					     function<map<int,double>(double const&, double const&)> const& InDistFunc,
					     double                                                  const& MuRef,
					     int                                                     const& PerturbativeOrder,
					     function<double(double const&)>                         const& Alphas,
					     int                                                     const& nsteps = 10);
}
