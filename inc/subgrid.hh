/*
  subgrid.hh:

  Author: Valerio Bertone
*/

#pragma once

using namespace std;

namespace apfel {

  /**
   * class subgrid:
   * Class for the x-space interpolation subgrids.
   * Subgrids are the building blocks of the interpolation procedure.
   * This class defines the "subgrid" object that includes, apart from the
   * grid itself, also the relevant parameters.
   **/
  class subgrid {
  private:
    // Attributes
    int    _nx;           // Number intervals
    double _xMin;         // Minumim value of x
    double _xMax;         // Maximum value of x (should always be 1)
    int    _InterDegree;  // Interpolation degree
    int    _xsize;        // Size of the subgrid
    bool   _IsExt;        // Is external
    double _Step;         // Step pf the logarthmically spaced grid
    double *_xsg;         // Actual grid

  public:
    // Constructors
    subgrid(int const& nx_, double const& xMin_, int const& InterDegree_); // Standard internal grid
    subgrid(int const& nx_, double *xsg_, int const& InterDegree_); // External grid
    subgrid(subgrid const& in_);

    // Getters
    int    nx()              const { return _nx; }
    int    InterDegree()     const { return _InterDegree; }
    double xMin()            const { return _xMin; }
    double xMax()            const { return _xMax; }
    bool   IsExt()           const { return _IsExt; }
    double Step()            const { return _Step; }
    double xg(int const& ix) const;

    // Interpolation functions
    double Interpolant(int const& beta, double const& x) const;

    // Check whether subgrids are equal 
    bool operator == (subgrid const& sg);
  };

}
