//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/constants.h"

namespace apfel
{
  /**
   * @brief Structure that contains all the precomputed quantities
   * needed to perform the DGLAP evolution in QCDxQED,
   * i.e. perturbative coefficients of splitting functions and
   * matching conditions, and the heavy-quark/lepton thresholds.
   */
  struct DglapObjectsQCDQED
  {
    double                                       Threshold;
    PartonSpecies                                Species;
    Set<Operator>                                UnitySet;
    std::map<std::pair<int, int>, Set<Operator>> SplittingFunctions;
    std::map<std::pair<int, int>, Set<Operator>> MatchingConditions;
  };

  /**
   * @name DGLAP object initializers in QCDxQED
   * Collection of functions that initialise DglapObjectsQED structure
   * for the different kinds of QCDxQED evolution currently available.
   */
  ///@{
  /**
   * @brief The InitializeDglapObjectsQCDQED function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions for QCDxQED evolution and store
   * them into a 'DglapObjectsQCDQED' structure.
   * @param g: the x-space grid
   * @param QuarkThresholds: the quark thresholds
   * @param LeptonThresholds: the lepton thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @param n3lo: whether N3LO corrections to splitting and matching functions are computer (default: false)
   * @param IMod: the vector of switches to vary the parameterisation of the approximated N3LO splitting functions (only relevant for n3lo = true) (default: all zero's)
   * @return A map of DglapObjectsQCDQED objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjectsQCDQED> InitializeDglapObjectsQCDQED(Grid                const& g,
                                                                 std::vector<double> const& QuarkThresholds,
                                                                 std::vector<double> const& LeptonThresholds,
                                                                 bool                const& OpEvol = false,
                                                                 double              const& IntEps = 1e-5,
                                                                 bool                const& n3lo = false,
                                                                 std::vector<int>    const& IMod = {0, 0, 0, 0, 0, 0, 0});
  ///@}

  /**
   * @name DGLAP builders for QCDxQED evolution
   * Collection of functions that build a Dglap object used to perform
   * the DGLAP QCDxQED evolution of distributions or operators.
   */
  ///@{
  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs the DGLAP QCDxQED evolution for distributions.
   * @param DglapObj: structure with the coefficients of the perturbative objects
   * @param InDistFunc: the distributions at the reference scale
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param Alphaem: the function returning the eletromagnetic coupling
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  std::unique_ptr<Dglap<Distribution>> BuildDglap(std::map<int, DglapObjectsQCDQED>                                  const& DglapObj,
                                                  std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc,
                                                  double                                                             const& MuRef,
                                                  int                                                                const& PerturbativeOrder,
                                                  std::function<double(double const&)>                               const& Alphas,
                                                  std::function<double(double const&)>                               const& Alphaem,
                                                  int                                                                const& nsteps = 10);

  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs the DGLAP QCDxQED evolution for operators.
   * @param DglapObj: structure with the coefficients of the perturbative objects
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param Alphaem: the function returning the eletromagnetic coupling
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  /*std::unique_ptr<Dglap<Operator>> BuildDglap(std::map<int, DglapObjectsQCDQED>    const& DglapObj,
                                              double                               const& MuRef,
                                              int                                  const& PerturbativeOrder,
                                              std::function<double(double const&)> const& Alphas,
                                              std::function<double(double const&)> const& Alphaem,
                                              int                                  const& nsteps = 10);*/
  ///@}
}
