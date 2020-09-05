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
    /// Flavour scheme
    enum FlavourScheme: int {VFNS, FFNS};

    /// Evolution theory
    enum EvolutionTheory: int {QCD, QCD_QED};

    /// Space- or time-like evolution (PDFs or FFs)
    enum Virtuality: int {SPACE, TIME};

    /// Polarisation of the evolution
    enum EvolPolarisation: int {UNP, POL, TRANS};

    /// Solution of the coupling RG equations
    enum CouplingEvolution: int {exact, expanded};

    /// Solution of the DGLAP equation
    enum PDFEvolution: int {exactmu, exactalpha, expandalpha, truncated};

    /// Heavy quark mass renormalisation
    enum MassRenScheme: int {POLE, MSBAR};
    ///@}

    /// Structure for the **subgrid** parameters
    struct GridParameters
    {
      int    nx;                  //!< Number of nodes
      double xmin;                //!< Minimum value of x
      int    id;                  //!< Interpolation degree
      std::vector<double> xgext;  //!< External grid (if any)
    };

    /**
     * @name Attributes of the struture.
     */
    ///@{
    std::string                 name;                //!< Identifier name of the setup
    double                      Q0;                  //!< Starting scale of the evolutions
    std::vector<GridParameters> GridParameters;      //!< Vector of the parameters of the subgrids
    bool                        Locked;              //!< Subgrid locked (ON, OFF)
    int                         nQg;                 //!< Number of the grid in Q
    double                      Qmin;                //!< Lower bound of the grid in Q
    double                      Qmax;                //!< Upper bound of the grid in Q
    int                         InterDegreeQ;        //!< Interpolation degree on the grid in Q
    double                      Lambda;              //!< Pole of the grid in ln(ln(Q/&Lambda;))
    int                         PerturbativeOrder;   //!< Perturbative order of the evolution (LO, NLO, NNLO)
    FlavourScheme               FlavourScheme;       //!< Flavour scheme (FFNS, VFNS)
    int                         Nf_FF;               //!< Number of active flavours in the FFNS
    EvolutionTheory             Theory;              //!< Evolution theory (QCD, QCD_QED)
    Virtuality                  Virtuality;          //!< Virtuality of the evolution (SPACE, TIME)
    EvolPolarisation            EvolPolarisation;    //!< Polarisation of the evolution (UNP, POL, TRANS)
    double                      AlphaQCDRef;         //!< Reference value of the QCD coupling
    double                      QQCDRef;             //!< reference scale of the QCD coupling
    double                      AlphaQEDRef;         //!< Reference value of the QED coupling
    double                      QQEDRef;             //!< Reference scale of the QCD coupling
    CouplingEvolution           CouplingEvolution;   //!< Solution of the RGE of the couplings
    PDFEvolution                PDFEvolution;        //!< Solution of the DGLAP equations
    double                      RenFacRatio;         //!< ratio muF / muR
    std::vector<double>         Thresholds;          //!< Heavy-quark thresholds
    std::vector<double>         Masses;              //!< Heavy-quark Masses
    MassRenScheme               MassRenScheme;       //!< Renormalization scheme for the heavy-quark masses (POLE, MSBAR)
    double                      TauMass;             //!< Mass of the &tau; lepton
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
      GridParameters{{100, 1e-5, 3, {}}, {60, 1e-1, 3, {}}, {50, 6e-1, 3, {}}, {50, 8e-1, 3, {}}},
      Locked(true),
      nQg(50), Qmin(1), Qmax(1000), InterDegreeQ(3), Lambda(0.25),
      PerturbativeOrder(2),
      FlavourScheme(VFNS),
      Nf_FF(3),
      Theory(QCD),
      Virtuality(SPACE),
      EvolPolarisation(UNP),
      AlphaQCDRef(0.35), QQCDRef(sqrt(2)),
      AlphaQEDRef(7.496252e-3), QQEDRef(1.777),
      CouplingEvolution(exact), PDFEvolution(exactmu),
      RenFacRatio(1),
      Thresholds{0, 0, 0, sqrt(2), 4.5, 175},
      Masses{Thresholds},
      MassRenScheme(POLE),
      TauMass(1.777),
      GaussAccuracy(1e-5),
      InSet({LHToyPDFs})
    {
    }
    // *INDENT-ON*
  };
}
