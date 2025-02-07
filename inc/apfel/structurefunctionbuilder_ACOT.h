/**
 * @file structurefunctionbuilderACOT.h
 * @author Peter Risse
 * @brief Same as structurefunctionbuilder.h but for the ACOT prescription
 * @version 1.0
 * @date 2024-07-31
 */

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"
#include "apfel/observable.h"
#include "apfel/disbasis.h"
#include "apfel/structurefunctionbuilder.h"

namespace apfel
{
  /**
   * @brief The InitializeF2NCObjectsSACOT precomputs the perturbative
   * coefficients of the coefficient functions for NC F2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF2NCObjectsSACOT(Grid const &g,
                                                                                                                  std::vector<double> const &Masses,
                                                                                                                  double const &IntEps = 1e-5,
                                                                                                                  int const &nQ = 150,
                                                                                                                  double const &Qmin = 1,
                                                                                                                  double const &Qmax = 300,
                                                                                                                  int const &intdeg = 3);
  /**
   * @brief The InitializeF2NCObjectsACOT precomputs the perturbative
   * coefficients of the coefficient functions for NC F2 in the full ACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF2NCObjectsACOT(Grid const &g,
                                                                                                                 std::vector<double> const &Masses,
                                                                                                                 double const &IntEps = 1e-5,
                                                                                                                 int const &nQ = 150,
                                                                                                                 double const &Qmin = 1,
                                                                                                                 double const &Qmax = 300,
                                                                                                                 int const &intdeg = 3);

  /**
   * @brief The InitializeF3NCObjectsSACOT precomputs the perturbative
   * coefficients of the coefficient functions for NC F3 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF3NCObjectsSACOT(Grid const &g,
                                                                                                                  std::vector<double> const &Masses,
                                                                                                                  double const &IntEps = 1e-5,
                                                                                                                  int const &nQ = 150,
                                                                                                                  double const &Qmin = 1,
                                                                                                                  double const &Qmax = 300,
                                                                                                                  int const &intdeg = 3);
  /**
   * @brief The InitializeFLNCObjectsSACOT precomputs the perturbative
   * coefficients of the coefficient functions for NC FL in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeFLNCObjectsSACOT(Grid const &g,
                                                                                                                  std::vector<double> const &Masses,
                                                                                                                  double const &IntEps = 1e-5,
                                                                                                                  int const &nQ = 150,
                                                                                                                  double const &Qmin = 1,
                                                                                                                  double const &Qmax = 300,
                                                                                                                  int const &intdeg = 3);
  /**
   * @brief The InitializeF2CCPlusObjectsSACOT precomputs the perturbative
   * coefficients of coefficient functions for combination ( F2(nu) +
   * F2(nubar) ) / 2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF2CCPlusObjectsSACOT(Grid const &g,
                                                                                                                      std::vector<double> const &Masses,
                                                                                                                      double const &IntEps = 1e-5,
                                                                                                                      int const &nQ = 150,
                                                                                                                      double const &Qmin = 1,
                                                                                                                      double const &Qmax = 300,
                                                                                                                      int const &intdeg = 3);
  /**
   * @brief The InitializeF2CCMinusObjectsSACOT precomputs the perturbative
   * coefficients of coefficient functions for combination ( F2(nu) -
   * F2(nubar) ) / 2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF2CCMinusObjectsSACOT(Grid const &g,
                                                                                                                       std::vector<double> const &Masses,
                                                                                                                       double const &IntEps = 1e-5,
                                                                                                                       int const &nQ = 150,
                                                                                                                       double const &Qmin = 1,
                                                                                                                       double const &Qmax = 300,
                                                                                                                       int const &intdeg = 3);
  /**
   * @brief The InitializeFLCCPlusObjectsSACOT precomputs the perturbative
   * coefficients of coefficient functions for combination ( FL(nu) +
   * FL(nubar) ) / 2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeFLCCPlusObjectsSACOT(Grid const &g,
                                                                                                                       std::vector<double> const &Masses,
                                                                                                                       double const &IntEps = 1e-5,
                                                                                                                       int const &nQ = 150,
                                                                                                                       double const &Qmin = 1,
                                                                                                                       double const &Qmax = 300,
                                                                                                                       int const &intdeg = 3);
  /**
   * @brief The InitializeFLCCMinusObjectsSACOT precomputs the perturbative
   * coefficients of coefficient functions for combination ( FL(nu) -
   * FL(nubar) ) / 2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeFLCCMinusObjectsSACOT(Grid const &g,
                                                                                                                       std::vector<double> const &Masses,
                                                                                                                       double const &IntEps = 1e-5,
                                                                                                                       int const &nQ = 150,
                                                                                                                       double const &Qmin = 1,
                                                                                                                       double const &Qmax = 300,
                                                                                                                       int const &intdeg = 3);
  /**
   * @brief The InitializeF3CCPlusObjectsSACOT precomputs the perturbative
   * coefficients of coefficient functions for combination ( F3(nu) +
   * F3(nubar) ) / 2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF3CCPlusObjectsSACOT(Grid const &g,
                                                                                                                       std::vector<double> const &Masses,
                                                                                                                       double const &IntEps = 1e-5,
                                                                                                                       int const &nQ = 150,
                                                                                                                       double const &Qmin = 1,
                                                                                                                       double const &Qmax = 300,
                                                                                                                       int const &intdeg = 3);
  /**
   * @brief The InitializeF3CCMinusObjectsSACOT precomputs the perturbative
   * coefficients of coefficient functions for combination ( F3(nu) -
   * F3(nubar) ) / 2 in the SACOT
   * scheme and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF3CCMinusObjectsSACOT(Grid const &g,
                                                                                                                       std::vector<double> const &Masses,
                                                                                                                       double const &IntEps = 1e-5,
                                                                                                                       int const &nQ = 150,
                                                                                                                       double const &Qmin = 1,
                                                                                                                       double const &Qmax = 300,
                                                                                                                       int const &intdeg = 3);
  /**
   * @brief The InitializeF2NCObjectsASACOT precomputs the perturbative
   * coefficients of coefficient functions for F2 in the aSACOT-chi(n) scheme up to
   * NNLO and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @param n: The scaling variable for the slow-rescaling with chi (default: 1)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF2NCObjectsASACOT(Grid const &g,
                                                                                                                   std::vector<double> const &Masses,
                                                                                                                   double const &IntEps = 1e-5,
                                                                                                                   int const &nQ = 100,
                                                                                                                   double const &Qmin = 1,
                                                                                                                   double const &Qmax = 300,
                                                                                                                   int const &intdeg = 3,
                                                                                                                   double const &n=1); 
  /**
   * @brief The InitializeFLNCObjectsASACOT precomputs the perturbative
   * coefficients of coefficient functions for FL in the aSACOT-chi(n) scheme up to
   * NNLO and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @param n: The scaling variable for the slow-rescaling with chi (default: 1)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeFLNCObjectsASACOT(Grid const &g,
                                                                                                                   std::vector<double> const &Masses,
                                                                                                                   double const &IntEps = 1e-5,
                                                                                                                   int const &nQ = 100,
                                                                                                                   double const &Qmin = 1,
                                                                                                                   double const &Qmax = 300,
                                                                                                                   int const &intdeg = 3,
                                                                                                                   double const &n=1); 

  /**
   * @brief The InitializeF3NCObjectsASACOT precomputs the perturbative
   * coefficients of coefficient functions for F3 in the aSACOT-chi(n) scheme up to
   * NNLO and stores them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Masses: the heavy quark masses
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @param nQ: the number of nodes of the grid in Q (default: 150)
   * @param Qmin: the lower bound of the grid in Q (default: 1)
   * @param Qmax: the upper bound of the grid in Q (default: 300)
   * @param intdeg: the interpolation degree on the grid in Q (default: 3)
   * @param n: The scaling variable for the slow-rescaling with chi (default: 1)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const &, std::vector<double> const &)> InitializeF3NCObjectsASACOT(Grid const &g,
                                                                                                                   std::vector<double> const &Masses,
                                                                                                                   double const &IntEps = 1e-5,
                                                                                                                   int const &nQ = 100,
                                                                                                                   double const &Qmin = 1,
                                                                                                                   double const &Qmax = 300,
                                                                                                                   int const &intdeg = 3,
                                                                                                                   double const &n=1); 
}