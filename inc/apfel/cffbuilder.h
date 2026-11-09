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
   * @brief The InitializeImCFF2NCObjectsZMSpaceLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * imaginary part of F2 for Q<SUP>2</SUP> &gt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFF2NCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeImCFF2NCObjectsZMTimeLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * imaginary part of F2 for Q<SUP>2</SUP> &lt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFF2NCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeReCFF2NCObjectsZMSpaceLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * real part of F2 for Q<SUP>2</SUP> &gt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFF2NCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeReCFF2NCObjectsZMTimeLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * real part of F2 for Q<SUP>2</SUP> &lt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFF2NCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeImCFFLNCObjectsZMSpaceLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * imaginary part of FL for Q<SUP>2</SUP> &gt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFFLNCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeImCFFLNCObjectsZMTimeLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * imaginary part of FL for Q<SUP>2</SUP> &lt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeImCFFLNCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeReCFFLNCObjectsZMSpaceLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * real part of FL for Q<SUP>2</SUP> &gt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFFLNCObjectsZMSpaceLike(Grid                const& g,
                                                                                                                          std::vector<double> const& Thresholds,
                                                                                                                          double              const& xi,
                                                                                                                          double              const& IntEps = 1e-5);

  /**
   * @brief The InitializeReCFFLNCObjectsZMTimeLike precomputes the
   * perturbative coefficients of coefficient functions for the NC
   * real part of FL for Q<SUP>2</SUP> &lt; 0 in the ZM scheme and stores
   * them in the 'StructureFunctionObjects' structure.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A StructureFunctionObjects-valued function
   */
  std::function<StructureFunctionObjects(double const&, std::vector<double> const&)> InitializeReCFFLNCObjectsZMTimeLike(Grid                const& g,
                                                                                                                         std::vector<double> const& Thresholds,
                                                                                                                         double              const& xi,
                                                                                                                         double              const& IntEps = 1e-5);
  ///@}
}
