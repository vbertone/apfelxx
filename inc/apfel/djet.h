//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @name Jet coefficients.
   * Perturbative coefficients for the definction of the low-scale jet
   * TMD depending on the jet algorithm.
   */
  ///@{
  /// &alpha;<SUB>s</SUB> correction to the cone-algorithm jet definition
  double dJetqCone1();

  /// &alpha;<SUB>s</SUB> correction to the kT-algorithm jet definition
  double dJetqkT1();
  ///@}
}
