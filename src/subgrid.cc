/*
 * APFEL++ 2017
 *
 * Authors: Valerio Bertone: valerio.bertone@vu.nl
 *          Stefano Carrazza: stefano.carrazza@cern.ch
 */

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

#include "apfel/subgrid.h"
#include "apfel/tools.h"

namespace apfel {

  //_________________________________________________________________________________
  SubGrid::SubGrid(int const& nx, double const& xMin, int const& InterDegree):
    _nx(nx),
    _xMin(xMin),
    _xMax(1),
    _InterDegree(InterDegree),
    _IsExternal(false)
  {
    // Compute grid
    _Step = log(_xMax / _xMin) / _nx;

    // building log spaced grid in x
    // number of points in x +1 (bins) + extra nodes for rhs interpolation.
    _xsg.resize(_nx+_InterDegree+1, 0);

    _xsg[0] = _xMin;
    const double exps = exp(_Step);
    for(auto ix = 1; ix < _xsg.size(); ix++) _xsg[ix] = _xsg[ix-1]*exps;
    _xsg[_nx] = 1;
  }

  //_________________________________________________________________________________
  SubGrid::SubGrid(vector<double> const& xsg, int const& InterDegree):
    _nx(xsg.size()-1),
    _xMin(xsg[0]),
    _xMax(1),
    _InterDegree(InterDegree),
    _IsExternal(true)
  {
    _xsg.resize(_nx+InterDegree+1, 0);
    copy(xsg.begin(), xsg.end(), _xsg.begin());

    // Check that the last point of the user-given grid is equal to one
    const double eps = 1e-12;
    if(fabs(_xsg[_nx]-1) >= eps)
      throw runtime_exception("SubGrid::SubGrid","The upper value of the external grid does not coincide with 1.");
    else
      _xsg[_nx] = 1;

    // Extend the grid for x > 1 for interpolation reasons using the same
    // width of the last bin in log scale
    const double step = log( xsg[_nx-1] / xsg[_nx-2] );
    const double exps = exp(step);
    for(auto ix = _nx; ix < _xsg.size(); ix++) _xsg[ix] = _xsg[ix-1] * exps;

    // In case on external grid the logarithmic step is not constant.
    // Set it to zero.
    _Step = 0;
  }

  //_________________________________________________________________________________
  double SubGrid::xg(int const& ix) const
  {
    assert(ix >= 0);
    assert(ix < _xsg.size());
    return _xsg[ix];
  }

  //_________________________________________________________________________________
  double SubGrid::Interpolant(int const& beta, double const& x) const
  {
    double w_int = 0;

    // Compute interpolant
    int bound = beta - _InterDegree;
    if(_InterDegree > beta) bound = 0;
    if(x < _xsg[bound] || x >= _xsg[beta+1]) return w_int;

    // If the grid is external ...
    if(_IsExternal) {
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

  //_________________________________________________________________________________
  bool SubGrid::operator == (SubGrid const& sg)
  {
    // Are both external or internal grids?
    if(_IsExternal != sg._IsExternal) return false;

    // In case they are external ...
    if(_IsExternal)
      {
        if(_nx != sg._nx)                   return false;
        if(_xsg.size() != sg._xsg.size())   return false;
        if(_xsg != sg._xsg)                 return false;
        if(_InterDegree != sg._InterDegree) return false;
      }
    else // In case they are internal ...
      {
        if(_xMin != sg._xMin)               return false;
        if(_xMax != sg._xMax)               return false;
        if(_InterDegree != sg._InterDegree) return false;
      }
    return true;
  }

  //_________________________________________________________________________________
  std::ostream& operator<<(std::ostream& os, const SubGrid& sg)
  {
    os << "SubGrid: " << &sg << "\n";
    os << "nx = " << sg._nx << "\n";
    os << "xMin = " << sg._xMin << "\n";
    os << "xMax = " << sg._xMax << "\n";
    os << "InterDegree = " << sg._InterDegree << "\n";
    os << "xsize = " << sg._xsg.size() << "\n";
    os << "IsExternal = " << sg._IsExternal << "\n";
    os << "Step = " << sg._Step << "\n";
    os << "xsg = ";
    for (const auto &v: sg._xsg)
      os << v << " ";
    os << "\n\n";
    return os;
  }
}
