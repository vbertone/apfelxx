//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include <iostream>
#include <algorithm>

#include "apfel/grid.h"
#include "apfel/subgrid.h"
#include "apfel/tools.h"

namespace apfel {

  //_________________________________________________________________________________
  bool ComparexMin(SubGrid const& sg1, SubGrid const& sg2)
  {
    if(sg1.xMin() == sg2.xMin())
      throw runtime_exception("ComparexMin","There are SubGrids with the same lower bound.");
    return sg1.xMin() < sg2.xMin();
  }

  //_________________________________________________________________________________
  Grid::Grid(vector<SubGrid> const& grs, bool const& lockgrids):
    _Locked(lockgrids),
    _ExtGrids(false),
    _GlobalGrid(grs)
  {
    CreateJointGrid();
  }

  //_________________________________________________________________________________
  void Grid::CreateJointGrid()
  {
    // Number of grids
    double const ng = _GlobalGrid.size();

    // Check if there are extenal grids and if so disable the locking
    for(int ig=0; ig<ng; ig++)
      if(_GlobalGrid[ig].IsExternal())
        _ExtGrids = true;

    if(_ExtGrids && _Locked)
      {
        warning("Grid::CreateJointedGrid", "External grids found... unlocking grids");
        _Locked = false;
      }

    // Now oder the SubGrids in such a way that they start with that having the lowest value of xMin
    // (only if there is more than one grid).
    if(ng > 1) sort(_GlobalGrid.begin(), _GlobalGrid.end(), ComparexMin);

    // In case there grids have been locked ...
    if(_Locked)
      {
        // Find the point of the "(ig-1)"-th SubGrid such that "x[ig-1][ix] < xMin[ig] < x[ig-1][ix+1]",
        // and replace "xMin[ig]" with "x[ig-1][ix]"
        for(auto ig = 1; ig < ng; ig++)
          {
            int const nxg     = _GlobalGrid[ig-1].nx();
            double const xmin = _GlobalGrid[ig].xMin();

            // Parameters of the adjusted grid
            int nx_new = -1;
            double xmin_new = -1;
            int const id_new = _GlobalGrid[ig].InterDegree();

            for(auto ix = 0; ix < nxg; ix++)
              {
                double const x = _GlobalGrid[ig-1].GetGrid()[ix];
                if(x > xmin)
                  {
                    xmin_new = x;
                    nx_new = nxg - ix;
                    break;
                  }
              }

            if (nx_new < 0 || xmin_new < 0)
              throw logic_exception("Grid::CreateJointedGrid", "SubGrids do not overlap.");

            // Find the closest multiple of "nx - ix + 1" to "nx",
            // i.e. "DensityFactor", and replace "nx" accordingly.
            int DensityFactor = _GlobalGrid[ig].nx() / nx_new;
            nx_new *= DensityFactor;

            // Compute the new SubGrid and replace it in the global grid
            SubGrid sgrid{nx_new, xmin_new, id_new};
            _GlobalGrid[ig] = sgrid;
          }
      }

    // Compute the joint grid.
    // Parameters of the joint grid
    int nx_joint = -1;
    int id_joint = _GlobalGrid[0].InterDegree(); // Use the interpolation degree of the first grid
    vector<double> xg_joint_vect;

    for(auto ig = 0; ig < ng; ig++)
      {
        int const nxg = _GlobalGrid[ig].nx();
        double xtrans;
        if(ig < ng-1) xtrans = _GlobalGrid[ig+1].xMin();
        else          xtrans = 1 + 2 * eps12;
        for(auto ix = 0; ix <= nxg; ix++)
          {
            double const x = _GlobalGrid[ig].GetGrid()[ix];
            if(xtrans - x < eps12) break;
            nx_joint++;
            xg_joint_vect.push_back(x);
          }
      }

    // initialize another SubGrid for the jointed grid.
    _JointGrid = unique_ptr<SubGrid>(new SubGrid{xg_joint_vect, id_joint});
  }

  //_________________________________________________________________________________
  bool Grid::operator == (Grid const& g) const
  {
    if(_Locked != g._Locked)                       return false;
    if(_ExtGrids != g._ExtGrids)                   return false;
    if(_GlobalGrid.size() != g._GlobalGrid.size()) return false;

    for(auto ig = 0 ; ig < (int) _GlobalGrid.size(); ig++)
      if(_GlobalGrid[ig] != g._GlobalGrid[ig])
        return false;

    return true;
  }

  //_________________________________________________________________________________
  std::ostream& operator<<(std::ostream& os, Grid const& gr)
  {
    os << "Grid: " << &gr << "\n";
    os << "Locked    = " << gr._Locked << "\n";
    os << "ExtGrids  = " << gr._ExtGrids << "\n";
    os << "JointGrid = " << &gr._JointGrid << "\n";
    for (const auto &v: gr._JointGrid->GetGrid()) os << v << " ";
    os << "\n\n";
    return os;
  }

}
