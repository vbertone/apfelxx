//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/grid.h"
#include "apfel/constants.h"
#include "apfel/messages.h"

namespace apfel
{
  //_________________________________________________________________________________
  Grid::Grid():
    _Locked(true),
    _ExtGrids(false),
    _GlobalGrid({})
  {
  }

  //_________________________________________________________________________________
  Grid::Grid(std::vector<SubGrid> const& grs, bool const& lockgrids):
    _Locked(lockgrids),
    _ExtGrids(false),
    _GlobalGrid(grs)
  {
    _JointGrid = CreateJointGrid();
  }

  //_________________________________________________________________________________
  bool ComparexMin(SubGrid const& sg1, SubGrid const& sg2)
  {
    if (sg1.xMin() == sg2.xMin())
      throw std::runtime_error(error("ComparexMin", "There are SubGrids with the same lower bound."));
    return sg1.xMin() < sg2.xMin();
  }

  //_________________________________________________________________________________
  SubGrid Grid::CreateJointGrid()
  {
    // Number of grids.
    double const ng = _GlobalGrid.size();

    // Check if there are extenal grids and if so disable the locking.
    for (int ig=0; ig<ng; ig++)
      if (_GlobalGrid[ig].IsExternal())
        _ExtGrids = true;

    if (_ExtGrids && _Locked)
      {
        warning("Grid::CreateJointGrid", "External grids found... unlocking grids");
        _Locked = false;
      }

    // Now oder the SubGrids in such a way that they start with that
    // having the lowest value of xMin (only if there is more than one
    // grid).
    if (ng > 1)
      std::sort(_GlobalGrid.begin(), _GlobalGrid.end(), ComparexMin);

    // In case there grids have been locked ...
    if (_Locked)
      // Find the point of the "(ig-1)"-th SubGrid such that
      // "x[ig-1][ix] < xMin[ig] < x[ig-1][ix+1]", and replace
      // "xMin[ig]" with "x[ig-1][ix]".
      for (int ig = 1; ig < ng; ig++)
        {
          const int nxg     = _GlobalGrid[ig-1].nx();
          const double xmin = _GlobalGrid[ig].xMin();

          // Parameters of the adjusted grid.
          int nx_new = -1;
          double xmin_new = -1;
          const int id_new = _GlobalGrid[ig].InterDegree();
          const std::vector<double> xg = _GlobalGrid[ig-1].GetGrid();

          for (int ix = 0; ix < nxg; ix++)
            if (xg[ix] > xmin)
              {
                xmin_new = xg[ix];
                nx_new = nxg - ix;
                break;
              }

          if (nx_new < 0 || xmin_new < 0)
            throw std::runtime_error(error("Grid::CreateJointGrid", "SubGrids do not overlap."));

          // Find the closest multiple of "nx - ix + 1" to "nx",
          // i.e. "DensityFactor", and replace "nx" accordingly.
          const int DensityFactor = _GlobalGrid[ig].nx() / nx_new;
          nx_new *= DensityFactor;

          // Compute the new SubGrid and replace it in the global
          // grid.
          _GlobalGrid[ig] = SubGrid{nx_new, xmin_new, id_new};
        }

    // Compute the joint grid. Parameters of the joint grid. Use the
    // interpolation degree of the first grid.
    const int id_joint = _GlobalGrid[0].InterDegree();
    std::vector<double> xg_joint_vect;

    for (int ig = 0; ig < ng; ig++)
      {
        const int nxg = _GlobalGrid[ig].nx();
        const std::vector<double> xg = _GlobalGrid[ig].GetGrid();

        double xtrans;
        if (ig < ng - 1)
          xtrans = _GlobalGrid[ig+1].xMin();
        else
          xtrans = 1 + 2 * eps12;
        for (int ix = 0; ix <= nxg; ix++)
          {
            if (xtrans - xg[ix] < eps12)
              break;
            xg_joint_vect.push_back(xg[ix]);
          }
      }

    // Initialize another SubGrid for the joint grid.
    return SubGrid{xg_joint_vect, id_joint};
  }

  //_________________________________________________________________________________
  bool Grid::operator == (Grid const& g) const
  {
    if (_Locked != g._Locked)
      return false;
    if (_ExtGrids != g._ExtGrids)
      return false;
    if (_GlobalGrid.size() != g._GlobalGrid.size())
      return false;

    for (int ig = 0 ; ig < (int) _GlobalGrid.size(); ig++)
      if (_GlobalGrid[ig] != g._GlobalGrid[ig])
        return false;

    return true;
  }

  //_________________________________________________________________________________
  bool Grid::operator != (Grid const& g) const
  {
    if (*this == g)
      return false;
    else
      return true;
  }

  //_________________________________________________________________________________
  std::ostream& operator<<(std::ostream& os, Grid const& gr)
  {
    os << "Grid: " << &gr << "\n";
    os << "Locked    = " << gr._Locked << "\n";
    os << "ExtGrids  = " << gr._ExtGrids << "\n";
    os << "JointGrid = " << &gr._JointGrid << "\n";
    for (const auto &v: gr._JointGrid.GetGrid()) os << v << " ";
    os << "\n\n";
    return os;
  }
}
