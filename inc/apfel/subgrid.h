/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone
 *          Stefano Carrazza
 */
#pragma once

namespace apfel
{
  /**
   * @brief Class for the x-space interpolation subgrids.
   *
   * Subgrids are the building blocks of the interpolation procedure.
   * This class defines the "subgrid" object that includes, apart from the
   * grid itself, also the relevant parameters.
   *
   */
  class subgrid
  {
  public:
    /**
     * @brief Standard subgrid constructor
     * @param nx_
     * @param xMin_
     * @param InterDegree_
     */
    subgrid(int const& nx_, double const& xMin_, int const& InterDegree_); // Standard internal grid

    /**
     * @brief External subgrid constructor
     * @param nx_
     * @param xsg_
     * @param InterDegree_
     */
    subgrid(int const& nx_, double *xsg_, int const& InterDegree_); // External grid

    /**
     * @brief Constructor to copy an existing subgrid
     * @param in_
     */
    subgrid(subgrid const& in_);

    // Getters
    int    nx()              const { return _nx; }
    int    InterDegree()     const { return _InterDegree; }
    bool   IsExt()           const { return _IsExt; }
    double xMin()            const { return _xMin; }
    double xMax()            const { return _xMax; }
    double Step()            const { return _Step; }

    /**
     * @brief Retrieve value of the subgrid in the "ix"-th point
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
     * @brief Check whether subgrids are equal
     * @param sg
     * @return
     */
    bool operator == (subgrid const& sg);

  private:
    int    _nx;           // Number intervals
    double _xMin;         // Minumim value of x
    double _xMax;         // Maximum value of x (should always be 1)
    int    _InterDegree;  // Interpolation degree
    int    _xsize;        // Size of the subgrid
    bool   _IsExt;        // Is external
    double _Step;         // Step pf the logarthmically spaced grid
    double *_xsg;         // Actual grid
  };
}
