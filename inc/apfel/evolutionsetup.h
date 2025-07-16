//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/lhtoypdfs.h"

#include <string>
#include <cmath>
#include <vector>
#include <functional>
#include <map>

namespace apfel
{
  /**
   * @brief The EvolutionSetup structure is a collection of all
   * possible evolution parameters.
   */
  struct EvolutionSetup
  {
    /**
     * @name Enumerators of the evolution setup structure
     */
    ///@{
    /// Evolution theory
    enum EvolutionTheory: int {QCD, QCD_QED};

    /// Space- or time-like evolution (PDFs or FFs)
    enum Virtuality: int {SPACE, TIME};

    /// Polarisation of the evolution
    enum EvolPolarisation: int {UNP, POL, TRANS};
    ///@}

    /// Structure for the **subgrid** parameters
    struct GridParameters
    {
      int    nx;                  //!< Number of nodes
      double xmin;                //!< Minimum value of x
      int    id;                  //!< Interpolation degree
    };

    /**
     * @name Attributes of the struture.
     */
    ///@{
    std::string                 name;                //!< Identifier name of the setup
    double                      Q0;                  //!< Starting scale of the evolutions
    std::vector<GridParameters> GridParameters;      //!< Vector of the parameters of the subgrids
    int                         nQg;                 //!< Number of the grid in Q
    double                      Qmin;                //!< Lower bound of the grid in Q
    double                      Qmax;                //!< Upper bound of the grid in Q
    int                         InterDegreeQ;        //!< Interpolation degree on the grid in Q
    int                         PerturbativeOrder;   //!< Perturbative order of the evolution (LO, NLO, NNLO)
    EvolutionTheory             Theory;              //!< Evolution theory (QCD, QCD_QED)
    Virtuality                  Virtuality;          //!< Virtuality of the evolution (SPACE, TIME)
    EvolPolarisation            EvolPolarisation;    //!< Polarisation of the evolution (UNP, POL, TRANS)
    double                      AlphaQCDRef;         //!< Reference value of the QCD coupling
    double                      AlphaQEDRef;         //!< Reference value of the QED coupling
    double                      QRef;                //!< Reference scale of the couplings
    std::vector<double>         QuarkThresholds;     //!< Heavy-quark thresholds
    std::vector<double>         QuarkMasses;         //!< Heavy-quark masses
    std::vector<double>         LeptonThresholds;    //!< Charged-lepton thresholds
    double                      GaussAccuracy;       //!< Accuracy of the dguass integration
    std::vector<std::function<std::map<int, double>(double const&, double const&)>> InSet; //!< Input set of distributions at the initial scale
    ///@}

    /**
     * @brief the constructor
     */
    // *INDENT-OFF*
    EvolutionSetup():
      name("default"),
      Q0(sqrt(2)),
      GridParameters{{100, 1e-5, 3}, {100, 1e-1, 3}, {100, 6e-1, 3}, {80, 8.5e-1, 5}},
      nQg(50), Qmin(1), Qmax(1000), InterDegreeQ(3),
      PerturbativeOrder(2),
      Theory(QCD),
      Virtuality(SPACE),
      EvolPolarisation(UNP),
      AlphaQCDRef(0.35),
      AlphaQEDRef(7.496252e-3),
      QRef(sqrt(2)),
      QuarkThresholds{0, 0, 0, sqrt(2), 4.5, 175},
      QuarkMasses{QuarkThresholds},
      LeptonThresholds{0, 0, 1.777},
      GaussAccuracy(1e-5),
      InSet({LHToyPDFs})
    {
    }
    // *INDENT-ON*
  };
}
