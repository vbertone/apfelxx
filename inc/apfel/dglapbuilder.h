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
    unordered_map<int,EvolutionBasisQCD> evbasis;
    unordered_map<int,MatchingBasisQCD>  matchbasis;
    unordered_map<int,Set<Operator>>     P0;
    unordered_map<int,Set<Operator>>     P1;
    unordered_map<int,Set<Operator>>     P2;
    unordered_map<int,Set<Operator>>     M0;
    unordered_map<int,Set<Operator>>     M1;
    unordered_map<int,Set<Operator>>     M2;
  };

  /**
   * @brief The InitializeDglapObjectsQCD, precompute the perturbative
   * coefficients of splitting functions and matching conditions and
   * store them in the 'DglapObjects' structure.
   * @param g the grid
   * @param IntEps the integration accuracy
   * @return
   */
  DglapObjects InitializeDglapObjectsQCD(Grid const& g, double const& IntEps = 1e-5);

  /**
   * @brief The DglapBuild, builds the dglap object
   *
   * @param DglapObj structure with the coefs. of perturbative objects
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param Masses the masses
   * @param Thresholds the threshold
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param nsteps the number of steps for RK.
   * @return
   */
  unique_ptr<Dglap> DglapBuild(DglapObjects                                                      const& DglapObj,
			       function<unordered_map<int,double>(double const&, double const&)> const& InDistFunc,
			       double                                                            const& MuRef,
			       vector<double>                                                    const& Masses,
			       vector<double>                                                    const& Thresholds,
			       int                                                               const& PerturbativeOrder,
			       function<double(double const&)>                                   const& Alphas,
			       int                                                               const& nsteps = 10);

  /**
   * @brief The DglapBuild, builds the dglap object
   *
   * @param DglapObj structure with the coefs. of perturbative objects
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param Masses the masses
   * @param Thresholds the threshold
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param nsteps the number of steps for RK.
   * @return
   */
  unique_ptr<Dglap> DglapBuild(DglapObjects                                                      const& DglapObj,
			       function<unordered_map<int,double>(double const&, double const&)> const& InDistFunc,
			       double                                                            const& MuRef,
			       vector<double>                                                    const& Masses,
			       int                                                               const& PerturbativeOrder,
			       function<double(double const&)>                                   const& Alphas,
			       int                                                               const& nsteps = 10);


  /**
   * @brief The DglapBuild, builds the dglap object
   *
   * @param DglapObj structure with the coefs. of perturbative objects
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param Masses the masses
   * @param Thresholds the threshold
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param nsteps the number of steps for RK.
   * @return
   */
  unique_ptr<Dglap> DglapBuild(DglapObjects                                               const& DglapObj,
			       function<double(int const&, double const&, double const&)> const& InDistFunc,
			       double                                                     const& MuRef,
			       vector<double>                                             const& Masses,
			       vector<double>                                             const& Thresholds,
			       int                                                        const& PerturbativeOrder,
			       function<double(double const&)>                            const& Alphas,
			       int                                                        const& nsteps = 10);

  /**
   * @brief The DglapBuild, builds the dglap object
   *
   * @param DglapObj structure with the coefs. of perturbative objects
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param Masses the masses
   * @param Thresholds the threshold
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param nsteps the number of steps for RK.
   * @return
   */
  unique_ptr<Dglap> DglapBuild(DglapObjects                                               const& DglapObj,
			       function<double(int const&, double const&, double const&)> const& InDistFunc,
			       double                                                     const& MuRef,
			       vector<double>                                             const& Masses,
			       int                                                        const& PerturbativeOrder,
			       function<double(double const&)>                            const& Alphas,
			       int                                                        const& nsteps = 10);
}
