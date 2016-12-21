/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@vu.nl
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */
#pragma once

#include <vector>
using std::vector;

namespace apfel
{
  /**
   * @brief Class for the x-space interpolation SubGrids.
   *
   * Subgrids are the building blocks of the interpolation procedure.
   * This class defines the "SubGrid" object that includes, apart from the
   * grid itself, also the relevant parameters.
   *
   */
  class SubGrid
  {
  public:

    SubGrid() = delete;

    /**
     * @brief Standard internal grid constructor.
     * @param nx number of grid points in x.
     * @param xMin lower edge x of the grid.
     * @param InterDegree interpolation degree.
     */
    SubGrid(int const& nx, double const& xMin, int const& InterDegree);

    /**
     * @brief External grid constructor
     * @param xsg a vector with the nodes of the grid
     * @param InterDegree interpolation degree
     */
    SubGrid(vector<double> const& xsg, int const& InterDegree);

    /**
     * @brief Retrieve value of the SubGrid in the "ix"-th point
     * @param ix
     * @return
     */
    double xg(int const& ix) const;

    /**
     * @brief Interpolant
     * @param beta
     * @param x
     * @return
     */
    double Interpolant(int const& beta, double const& x) const;

    /**
     * @brief Check whether SubGrids are equal
     * @param sg
     * @return
     */
    bool operator == (SubGrid const& sg);

    // Getters
    int    nx()          const { return _nx; }          //!< return the number of x points
    int    InterDegree() const { return _InterDegree; } //!<
    bool   IsExternal()  const { return _IsExternal; }  //!<
    double xMin()        const { return _xMin; }        //!<
    double xMax()        const { return _xMax; }        //!<
    double Step()        const { return _Step; }        //!< return the step size of the log grid

  private:
    int    _nx;           // Number intervals
    double _xMin;         // Minumim value of x
    double _xMax;         // Maximum value of x (should always be 1)
    int    _InterDegree;  // Interpolation degree
    bool   _IsExternal;   // Is external
    double _Step;         // Step pf the logarthmically spaced grid
    vector<double> _xsg;  // Actual grid

    friend std::ostream& operator<<(std::ostream& os, const SubGrid& dt);
  };

  /**
   * @brief operator <<
   * @param os
   * @param dt
   * @return
   */
  std::ostream& operator<<(std::ostream& os, const SubGrid& dt);
}
