//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/dglapnonlinear.h"

#include <memory>

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
    std::map<int, Set<Operator>>     SplittingFunctions;
    std::map<int, Set<Operator>>     MatchingConditions;
    std::map<int, Set<Distribution>> InhomogeneousTerms;
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
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @param IMod: the vector of switches to vary the parameterisation of the approximated N3LO splitting functions (default: all zero)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCD(Grid                const& g,
                                                        std::vector<double> const& Masses,
                                                        std::vector<double> const& Thresholds,
                                                        bool                const& OpEvol = false,
                                                        double              const& IntEps = 1e-5,
                                                        std::vector<int>    const& IMod = {0, 0, 0, 0, 0, 0, 0});

  /**
   * @brief The InitializeDglapObjectsQCD function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @param IMod: the vector of switches to vary the parameterisation of the approximated N3LO splitting functions (default: all zero)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCD(Grid                const& g,
                                                        std::vector<double> const& Thresholds,
                                                        bool                const& OpEvol = false,
                                                        double              const& IntEps = 1e-5,
                                                        std::vector<int>    const& IMod = {0, 0, 0, 0, 0, 0, 0});

  /**
   * @brief The InitializeDglapObjectsQCD function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions with MSbar masses and store
   * them into a 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDMSbarMass(Grid                const& g,
                                                                 std::vector<double> const& Masses,
                                                                 std::vector<double> const& Thresholds,
                                                                 bool                const& OpEvol = false,
                                                                 double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCD function precomputes the
   * perturbative coefficients of space-like unpolarised splitting
   * functions and matching conditions with MSbar masses and store
   * them into a 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark masses
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDMSbarMass(Grid                const& g,
                                                                 std::vector<double> const& Thresholds,
                                                                 bool                const& OpEvol = false,
                                                                 double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDpol function precomputes the
   * perturbative coefficients of space-like longitudinally polarised
   * splitting functions and matching conditions (assumed to be equal
   * to the unpolarised ones) and store them into a 'DglapObjects'
   * structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDpol(Grid                const& g,
                                                           std::vector<double> const& Masses,
                                                           std::vector<double> const& Thresholds,
                                                           bool                const& OpEvol = false,
                                                           double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDpol function precomputes the
   * perturbative coefficients of space-like longitudinally polarised
   * splitting functions and matching conditions (assumed to be equal
   * to the unpolarised ones) and store them into a 'DglapObjects'
   * structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDpol(Grid                const& g,
                                                           std::vector<double> const& Thresholds,
                                                           bool                const& OpEvol = false,
                                                           double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDT function precomputes the
   * perturbative coefficients of time-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDT(Grid                const& g,
                                                         std::vector<double> const& Masses,
                                                         std::vector<double> const& Thresholds,
                                                         bool                const& OpEvol = false,
                                                         double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDT function precomputes the
   * perturbative coefficients of time-like unpolarised splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDT(Grid                const& g,
                                                         std::vector<double> const& Thresholds,
                                                         bool                const& OpEvol = false,
                                                         double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDTpol function precomputes the
   * perturbative coefficients of time-like longitudinally polarised
   * splitting functions and matching conditions and store them into a
   * 'DglapObjects' structure. Limited to LO for now.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDTpol(Grid                const& g,
                                                            std::vector<double> const& Masses,
                                                            std::vector<double> const& Thresholds,
                                                            bool                const& OpEvol = false,
                                                            double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDTpol function precomputes the
   * perturbative coefficients of time-like longitudinally polarised
   * splitting functions and matching conditions and store them into a
   * 'DglapObjects' structure. Limited to LO for now.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDTpol(Grid                const& g,
                                                            std::vector<double> const& Thresholds,
                                                            bool                const& OpEvol = false,
                                                            double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of space-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDtrans(Grid                const& g,
                                                             std::vector<double> const& Masses,
                                                             std::vector<double> const& Thresholds,
                                                             bool                const& OpEvol = false,
                                                             double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of space-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDtrans(Grid                const& g,
                                                             std::vector<double> const& Thresholds,
                                                             bool                const& OpEvol = false,
                                                             double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of timelike-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy-quark masses
   * @param Thresholds: the heavy quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDTtrans(Grid                const& g,
                                                              std::vector<double> const& Masses,
                                                              std::vector<double> const& Thresholds,
                                                              bool                const& OpEvol = false,
                                                              double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeDglapObjectsQCDtrans function precomputes
   * the perturbative coefficients of time-like transverity splitting
   * functions and matching conditions and store them into a
   * 'DglapObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy-quark thresholds
   * @param OpEvol: the switch for the computation of the evolution operator (default: false)
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   * @note This function assumes that masses and thresholds coincide.
   */
  std::map<int, DglapObjects> InitializeDglapObjectsQCDTtrans(Grid                const& g,
                                                              std::vector<double> const& Thresholds,
                                                              bool                const& OpEvol = false,
                                                              double              const& IntEps = 1e-5);
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
   * @param xi: the scale-variation parameter (default: 1)
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  std::unique_ptr<Dglap<Distribution>> BuildDglap(std::map<int, DglapObjects>                                        const& DglapObj,
                                                  std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc,
                                                  double                                                             const& MuRef,
                                                  int                                                                const& PerturbativeOrder,
                                                  std::function<double(double const&)>                               const& Alphas,
                                                  double                                                             const& xi = 1,
                                                  int                                                                const& nsteps = 10);

  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs the DGLAP evolution for operators.
   * @param DglapObj: structure with the coefficients of the perturbative objects
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param xi: the scale-variation parameter (default: 1)
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  std::unique_ptr<Dglap<Operator>> BuildDglap(std::map<int, DglapObjects>          const& DglapObj,
                                              double                               const& MuRef,
                                              int                                  const& PerturbativeOrder,
                                              std::function<double(double const&)> const& Alphas,
                                              double                               const& xi = 1,
                                              int                                  const& nsteps = 10);

  /**
    * @brief The BuildDglap function builds the actual dglap object
    * that performs the DGLAP evolution for distributions.
    * @param DglapObj: DglapObjects-valued function that returns the structure with the coefficients of the perturbative objects as functions of the scale
    * @param Thresholds: the heavy-quark thresholds
    * @param InDistFunc: the distributions at the reference scale
    * @param MuRef: the reference scale
    * @param PerturbativeOrder: the perturbative order of the evolution
    * @param Alphas: the function returning the strong coupling
    * @param nsteps: the number of steps of the ODE solver (default: 10).
    * @return A unique pointer to a Dglap object
    */
  std::unique_ptr<Dglap<Distribution>> BuildDglap(std::function<DglapObjects(double const&)>                         const& DglapObj,
                                                  std::vector<double>                                                const& Thresholds,
                                                  std::function<std::map<int, double>(double const&, double const&)> const& InDistFunc,
                                                  double                                                             const& MuRef,
                                                  int                                                                const& PerturbativeOrder,
                                                  std::function<double(double const&)>                               const& Alphas,
                                                  int                                                                const& nsteps = 10);

  /**
   * @brief The BuildDglap function builds the actual dglap object
   * that performs a non-liner DGLAP evolution.
   * @param DglapObj: structure with the coefficients of the perturbative objects
   * @param TranformationFuncs: set of functions that trasform distributions in the l.h.s. of DGLAP
   * @param InDistFunc: the distributions at the reference scale
   * @param MuRef: the reference scale
   * @param PerturbativeOrder: the perturbative order of the evolution
   * @param Alphas: the function returning the strong coupling
   * @param xi: the scale-variation parameter (default: 1)
   * @param nsteps: the number of steps of the ODE solver (default: 10).
   * @return A unique pointer to a Dglap object
   */
  std::unique_ptr<DglapNonLinear> BuildDglapNonLinear(std::map<int, DglapObjects>                                                                      const& DglapObj,
                                                      std::function<std::map<int, std::function<double(std::map<int, double> const&)>>(double const&)> const& TranformationFuncs,
                                                      std::function<std::map<int, double>(double const&, double const&)>                               const& InDistFunc,
                                                      double                                                                                           const& MuRef,
                                                      int                                                                                              const& PerturbativeOrder,
                                                      std::function<double(double const&)>                                                             const& Alphas,
                                                      double                                                                                           const& xi = 1,
                                                      int                                                                                              const& nsteps = 10);
  ///@}
}
