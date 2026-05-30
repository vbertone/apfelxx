//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/dglapbuilder.h"

namespace apfel
{

  /**
   * @name Factorisation scheme changers
   * Suite of functions that allow to change the factoristion scheme
   * of a set of DglapObjects given a set of kernels.
   */
  ///@{
  /**
   * @brief The InitializeSchemChangeKernelsKrk function precomputes
   * the perturbative coefficients of of the Krk to MSbar
   * scheme-change kernels and store them into a A map of maps of
   * Set<Operator> objects.
   * @param g: the x-space grid
   * @param Thresholds: the quark thresholds
   * @param IntEps: the integration accuracy (default: 10<SUP>-5</SUP>)
   * @return A map of maps of Set<Operator> objects, one for each
   * possible nf and perturbative order
   */
  std::map<int, std::map<int, Set<Operator>>> InitializeSchemeChangeKernelsKrk(Grid                const& g,
                                                                               std::vector<double> const& Thresholds,
                                                                               double              const& IntEps = 1e-5);

  /**
   * @brief The ChangeFactorisationSchemeMSbarToK function changes the
   * factoristion scheme from MSbar to the scheme coded in the input
   * kernels K.
   * @param DOMSbar: DGLAP objects in the MSbar scheme
   * @param K: scheme-change kernels
   */
  std::map<int, DglapObjects> ChangeFactorisationSchemeMSbarToK(std::map<int, DglapObjects>                 const& DOMSbar,
                                                                std::map<int, std::map<int, Set<Operator>>> const& K);
  ///@}
}
