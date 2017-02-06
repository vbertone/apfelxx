//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <array>
#include <apfel/qgrid.h>

using std::array;

namespace apfel
{
  /**
   * @brief Alpha_s on a QGrid:
   *
   * This class fills out the _GridValues vector step by step as the
   * calculation of the next step is always computed by starting from
   * the previous step. This method is performing because the steps
   * are typically small and thus it is sufficient to make one single
   * step in the Runge-Kutta algorithm.
   */
  class GridAlphaQCD: public QGrid<double>
  {
  public:

    GridAlphaQCD() = delete;

    /**
     * @brief GridAlphaQCD default constructor.
     * @param AlphaRef the reference value of the coupling.
     * @param MuRef the reference value of the scale.
     * @param Masses vector of masses.
     * @param Thresholds vector of thresholds.
     * @param pt perturbative order.
     * @param nQ number of Q nodes onbthe QGrid.
     * @param Qmin minimum value of the QGrid.
     * @param Qmax maximun value of the QGrid.
     * @param InterDegree interpolation degree.
     */
    GridAlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, vector<double> const& Thresholds, int const& pt, int const& nQ, double const& QMin, double const& QMax, int const& InterDegree);

    /**
     * @brief GridAlphaQCD constructor. Same as above but assuming Masses = Thresholds.
     */
    GridAlphaQCD(double const& AlphaRef, double const& MuRef, vector<double> const& Masses, int const& pt, int const& nQ, double const& QMin, double const& QMax, int const& InterDegree);

  };

}
