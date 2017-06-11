//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/dglap.h"

#include <functional>
#include <vector>

using std::function;
using std::vector;

namespace apfel
{
  /**
   * @brief The DglapBuildQCD, builds the dglap object for QCD
   *
   * @param g the grid
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param Masses the masses
   * @param Thresholds the threshold
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param IntEps the integration accuracy
   * @param nsteps the number of steps for RK.
   * @return
   */
  Dglap DglapBuildQCD(Grid                                                       const& g,
                      function<double(int const&, double const&, double const&)> const& InDistFunc,
                      double                                                     const& MuRef,
                      vector<double>                                             const& Masses,
                      vector<double>                                             const& Thresholds,
                      int                                                        const& PerturbativeOrder,
                      function<double(double const&)>                            const& Alphas,
                      double                                                     const& IntEps = 1e-5,
                      int                                                        const& nsteps = 10);

  /**
   * @brief The DglapBuildQCD, builds the dglap object for QCD
   *
   * @param g the grid
   * @param InDistFunc the PDF method to query flavors
   * @param MuRef the reference scale
   * @param Masses the masses
   * @param Thresholds the threshold
   * @param PerturbativeOrder the perturbative order
   * @param Alphas the alpha strong object
   * @param IntEps the integration accuracy
   * @param nsteps the number of steps for RK.
   * @return
   */
  Dglap DglapBuildQCD(Grid                                                       const& g,
                      function<double(int const&, double const&, double const&)> const& InDistFunc,
                      double                                                     const& MuRef,
                      vector<double>                                             const& Masses,
                      int                                                        const& PerturbativeOrder,
                      function<double(double const&)>                            const& Alphas,
                      double                                                     const& IntEps = 1e-5,
                      int                                                        const& nsteps = 10);

}
