/*
  subgrid.cc:

  Author: Valerio Bertone
*/

#include <iostream>
#include <cmath>
#include <stdexcept>

#include "../inc/subgrid.hh"

using namespace std;

namespace apfel {

  // ================================================================================
  // Standard subgrid constructor
  subgrid::subgrid(int const& nx_, double const& xMin_, int const& InterDegree_):
    _nx(nx_),
    _xMin(xMin_),
    _xMax(1),
    _InterDegree(InterDegree_),
    _xsize(_nx+_InterDegree+1),
    _IsExt(false)
  {
    // Compute grid
    _Step   = log(_xMax / _xMin) / _nx;
    _xsg    = new double[_xsize];
    _xsg[0] = _xMin;
    for(int ix=1; ix<_xsize; ix++) _xsg[ix] = _xsg[ix-1] * exp( _Step );

    // to avoid numerical problems, set _xsg[_nx] = 1 by hand
    _xsg[_nx] = 1;
  }

  // ================================================================================
  // External subgrid constructor
  subgrid::subgrid(int const& nx_, double *xsg_, int const& InterDegree_):
    _nx(nx_-1),
    _xMin(xsg_[0]),
    _xMax(1),
    _InterDegree(InterDegree_),
    _xsize(_nx+_InterDegree+1),
    _IsExt(true)
  {
    // Check that the last point of the user-given grid is equal to one
    double eps = 1e-12;
    if(abs(xsg_[_nx]-1) >= eps) throw runtime_error("The upper value of the external grid does not coincide with 1.");
    else                        xsg_[_nx] = 1;
      
    // Import grid
    _xsg = new double[_xsize];
    for(int ix=0; ix<_nx; ix++) _xsg[ix] = xsg_[ix];

    // Extend the grid for x > 1 for interpolation reasons using the same
    // width of the last bin in log scale
    double step = log( xsg_[_nx-1] / xsg_[_nx-2] );
    for(int ix=_nx; ix<_xsize; ix++) _xsg[ix] = _xsg[ix-1] * exp( step );

    // In case on external grid the logarithmic step is not constant.
    // Set it to zero.
    _Step = 0;
  }

  // ================================================================================
  // Constructor to copy an existing subgrid
  subgrid::subgrid(subgrid const& in_):
    _nx(in_._nx),
    _xMin(in_._xMin),
    _xMax(in_._xMax),
    _InterDegree(in_._InterDegree),
    _xsize(in_._xsize),
    _IsExt(in_._IsExt),
    _Step(in_._Step),
    _xsg(in_._xsg)
  { }

  // ================================================================================
  // Retrieve value of the subgrid in the "ix"-th point
  double subgrid::xg(int const& ix) const
  {
    if(ix < 0 || ix >= _xsize) throw runtime_error("Grid index out of range.");
    return _xsg[ix];
  }

  // ================================================================================
  // Interpolation functions
  double subgrid::Interpolant(int const& beta, double const& x) const
  {
    double w_int = 0;

    // Compute interpolant
    int bound = beta - _InterDegree;
    if(_InterDegree > beta) bound = 0;
    if(x < _xsg[bound] || x >= _xsg[beta+1]) return w_int;

    // If the grid is external ...
    if(_IsExt) {
      w_int = 1;
      int j;
      for(j=0; j<=beta-bound; j++) if(x >= _xsg[beta-j] && x < _xsg[beta-j+1]) break;
      for(int delta=0; delta<=_InterDegree; delta++) if(delta != j) w_int = w_int * log( x / _xsg[beta-j+delta] )
								      / log( _xsg[beta] / _xsg[beta-j+delta] );
    }
    // If the grid is internal ...
    else {
      w_int = 1;
      double fact = log( x / _xsg[beta] ) / _Step;
      int j;
      for(j=0; j<=beta-bound; j++) if(x >= _xsg[beta-j] && x < _xsg[beta-j+1]) break;
      for(int delta=0; delta<=_InterDegree; delta++) if(delta != j) w_int = w_int * ( fact / ( j - delta ) + 1 );
    }
    return w_int;
  }

  // ================================================================================
  // Check whether subgrids are equal
  bool subgrid::operator == (subgrid const& sg)
  {
    // Are both external or internal grids?                                                                                                                       
    if(_IsExt != sg._IsExt) return false;

    // In case they are external ...                                                                                                                              
    if(_IsExt) {
      if(_nx != sg._nx)                   return false;
      for(int ix=0; ix<_nx; ix++) {
	if(_xsg[ix] != sg._xsg[ix])       return false;
      }
      if(_InterDegree != sg._InterDegree) return false;
    }
    // In case they are internal ...                                                                                                                              
    else {

      if(_xMin != sg._xMin)               return false;
      if(_xMax != sg._xMax)               return false;
      if(_InterDegree != sg._InterDegree) return false;
    }
    return true;
  }

}
