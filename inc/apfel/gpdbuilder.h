//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/dglapbuilder.h"

namespace apfel
{
  /**
   * @name GPD object initializers
   * Collection of functions that initialise DglapObjects structure
   * for the different kinds of GPD evolution currently available.
   */
  ///@{
  /**
   * @brief The InitializeGPDObjects function precomputes the
   * perturbative coefficients of unpolarised GPD evolution kernels
   * and store them into a 'DglapObjects' structure. GPDs are assumed
   * to be continuous over heavy-quark thresholds.
   * @param g: the x-space grid
   * @param Thresholds: the heavy quark thresholds
   * @param xi: value of the skewness
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of DglapObject objects, one for each possible nf
   */
  std::map<int, DglapObjects> InitializeGpdObjects(Grid                const& g,
                                                   std::vector<double> const& Thresholds,
                                                   double              const& xi,
                                                   double              const& IntEps = 1e-5);
  ///@}
}
