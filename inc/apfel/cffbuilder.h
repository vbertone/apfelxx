//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/structurefunctionbuilder.h"

namespace apfel
{
  /**
   * @name Compton-form-factor object initializers
   * Collection of functions that initialise StructureFunctionObjects
   * structure for the different kinds of Compton form factors
   * available.
   */
  ///@{
  /**
   * @brief The InitializeImCFF1NCObjectsZM precomputes the
   * perturbative coefficients of coefficient functions for NC
   * imaginary part of F1 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.

   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFF1NCObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& xi,
                                                                                                                 double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeReCFF1NCObjectsZM precomputes the
   * perturbative coefficients of coefficient functions for NC real
   * part of F1 in the ZM scheme and store them in the
   * 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>})
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFF1NCObjectsZM(Grid                const& g,
                                                                                                                 std::vector<double> const& Thresholds,
                                                                                                                 double              const& xi,
                                                                                                                 double              const& IntEps = 1e-5);
  ///@}
}
