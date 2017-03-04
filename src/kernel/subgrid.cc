//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/subgrid.h"
#include "apfel/tools.h"

#include <iostream>
#include <cmath>
#include <cassert>
#include <algorithm>

namespace apfel
{
  //_________________________________________________________________________________
  SubGrid::SubGrid(int const& nx, double const& xMin, int const& InterDegree):
    _nx(nx),
    _InterDegree(InterDegree),
    _IsExternal(false),
    _xMin(xMin),
    _xMax(1)
  {
    // Compute grid
    _Step = log(_xMax / _xMin) / _nx;

    // building log spaced grid in x
    // number of points in x + 1 (bins) + extra nodes for rhs interpolation.
    _xsg.resize(_nx+_InterDegree+1, 0);

    _xsg[0] = _xMin;
    double const exps = exp(_Step);
    for(auto ix = 1; ix < (int) _xsg.size(); ix++) _xsg[ix] = _xsg[ix-1] * exps;
    _xsg[_nx] = 1;

    _lxsg.resize(_xsg.size());
    for (auto ix = 0; ix < (int) _xsg.size(); ix++) _lxsg[ix] = log(_xsg[ix]);
  }

  //_________________________________________________________________________________
  SubGrid::SubGrid(vector<double> const& xsg, int const& InterDegree):
    _nx(xsg.size()-1),
    _InterDegree(InterDegree),
    _IsExternal(true),
    _xMin(xsg[0]),
    _xMax(1),
    _Step(0)
  {
    _xsg.resize(_nx+InterDegree+1, 0);
    copy(xsg.begin(), xsg.end(), _xsg.begin());

    // Check that the last point of the user-given grid is equal to one
    if(fabs(_xsg[_nx]-1) >= eps11)
      throw runtime_exception("SubGrid::SubGrid","The upper value of the external grid does not coincide with 1.");
    else
      _xsg[_nx] = 1;

    // Extend the grid for x > 1 for interpolation reasons using the same
    // width of the last bin in log scale
    double const step = - log( xsg[_nx-1] );
    double const exps = exp(step);
    for(auto ix = _nx; ix < (int) _xsg.size(); ix++) _xsg[ix] = _xsg[ix-1] * exps;

    _lxsg.resize(_xsg.size());
    for (auto ix = 0; ix < (int) _xsg.size(); ix++) _lxsg[ix] = log(_xsg[ix]);
  }

  //_________________________________________________________________________________
  bool SubGrid::operator == (SubGrid const& sg) const
  {
    // Are both external or internal grids?
    if(_IsExternal != sg._IsExternal) return false;

    // In case they are external ...
    if(_IsExternal)
      {
        if(_xsg != sg._xsg)                 return false;
        if(_InterDegree != sg._InterDegree) return false;
      }
    else // In case they are internal ...
      {
        if(_nx != sg._nx)                   return false;
        if(_xMin != sg._xMin)               return false;
        if(_xMax != sg._xMax)               return false;
        if(_InterDegree != sg._InterDegree) return false;
      }
    return true;
  }

  //_________________________________________________________________________________
  bool SubGrid::operator != (SubGrid const& sg) const
  {
    if (*this == sg) return false;
    else             return true;
  }

  //_________________________________________________________________________________
  SubGrid& SubGrid::operator = (SubGrid const& sg)
  {
    // Copy attributes
    _IsExternal = sg._IsExternal;
    if(_IsExternal)
      {
        _xsg         = sg._xsg;
        _InterDegree = sg._InterDegree;
      }
    else
      {
        _nx          = sg._nx;
        _xMin        = sg._xMin;
        _xMax        = sg._xMax;
        _InterDegree = sg._InterDegree;
      }
    return *this;
  }

  //_________________________________________________________________________________
  std::ostream& operator<<(std::ostream& os, SubGrid const& sg)
  {
    os << "SubGrid: " << &sg << "\n";
    os << "nx          = " << sg._nx << "\n";
    os << "xMin        = " << sg._xMin << "\n";
    os << "xMax        = " << sg._xMax << "\n";
    os << "InterDegree = " << sg._InterDegree << "\n";
    os << "xsize       = " << sg._xsg.size() << "\n";
    os << "IsExternal  = " << sg._IsExternal << "\n";
    os << "Step        = " << sg._Step << "\n";
    os << "xsg         = ";
    for (const auto &v: sg._xsg) os << v << " ";
    os << "\n\n";
    return os;
  }
}
