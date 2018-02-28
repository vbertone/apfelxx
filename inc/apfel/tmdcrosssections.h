//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/set.h"
#include "apfel/distribution.h"

namespace apfel
{
  /**
   * @name TMD Drell-Yan cross sections
   * Functions that return the fully differential Drell-Yan cross
   * section as function of the transverse momentum on the vector
   * boson qT, taking as inputs the evolved TMD PDFs, the
   * center-of-mass energy Vs, the vector-boson virtuality Q and
   * rapidity, the final TMD scales muf and zetaf.
   */
  ///@{
  std::function<double(double const&)> TmdCrossSectionDY(double                                                                        const& Vs,
							 double                                                                        const& Q,
							 double                                                                        const& y,
							 std::function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDPDFs,
							 std::function<double(double const&)>                                          const& Alphas,
							 std::function<std::vector<double>(double const&)>                             const& fEWCharges,
							 int                                                                           const& PerturbativeOrder,
							 std::vector<double>                                                           const& Thresholds,
							 double                                                                        const& cmuf = 1,
							 double                                                                        const& czetaf = 1);

  //_____________________________________________________________________________
  std::function<double(double const&)> TmdCrossSectionDY(double                                                                          const& Vs,
							 double                                                                          const& Qmin,
							 double                                                                          const& Qmax,
							 double                                                                          const& ymin,
							 double                                                                          const& ymax,
							 std::function<Set<Distribution>(double const&)>                                 const& InTMDPDFs,
							 std::function<std::vector<double>(double const&, double const&, double const&)> const& EvolFact,
							 std::function<double(double const&)>                                            const& Alphas,
							 std::function<std::vector<double>(double const&)>                               const& fEWCharges,
							 int                                                                             const& PerturbativeOrder,
							 std::vector<double>                                                             const& Thresholds,
							 double                                                                          const& cmuf = 1,
							 double                                                                          const& czetaf = 1,
							 double                                                                          const& IntEps = 1e-5);
  ///@}

  /**
   * @name SIDIS TMD cross sections
   * Functions that return the fully differential SIDIS cross section
   * as function of the transverse momentum on the vector boson qT,
   * taking as inputs the evolved TMD PDFs and FFs, the center-of-mass
   * energy Vs, the inelasticity y, the momentum fractions of the
   * incmoming and outgoing hadrons x and z, rapidity, the final TMD
   * scales muf and zetaf.
   */
  ///@{
  std::function<double(double const&)> TmdCrossSectionSIDIS(double                                                                        const& Vs,
							    double                                                                        const& x,
							    double                                                                        const& y,
							    double                                                                        const& z,
							    std::function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDPDFs,
							    std::function<Set<Distribution>(double const&, double const&, double const&)> const& EvolvedTMDFFs,
							    std::function<double(double const&)>                                          const& Alphas,
							    std::function<std::vector<double>(double const&)>                             const& fEWCharges,
							    int                                                                           const& PerturbativeOrder,
							    std::vector<double>                                                           const& Thresholds,
							    double                                                                        const& cmuf = 1,
							    double                                                                        const& czetaf = 1);
  ///@}

  /**
   * @name Form factors
   */
  ///@{
 /**
   * @brief Perturbative form factor for Drell-Yan.
   */
  double HardCrossSectionDY(int const& pt, double const& Alphas, int const& nf, double const& kappa);

 /**
   * @brief Perturbative form factor for SIDIS.
   */
  double HardCrossSectionSIDIS(int const& pt, double const& Alphas, int const& nf, double const& kappa);
  ///@}
}
