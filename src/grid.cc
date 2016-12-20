/*
  grid.cc:

  Author: Valerio Bertone
*/

#include <iostream>
#include <stdexcept>
#include <algorithm>

#include "apfel/grid.h"
#include "apfel/subgrid.h"

using namespace std;

namespace apfel {

  // ================================================================================
  // Constructor
  grid::grid():
    _Locked(false),
    _ExtGrids(false)
  { }

  // ================================================================================
  // Destructor
  grid::~grid()
  {
    _GlobalGrid.erase(_GlobalGrid.begin(), _GlobalGrid.end());
  }

  // ================================================================================
  // Add a SubGrid
  void grid::AddGrid(SubGrid const& sgrid_) { _GlobalGrid.push_back(sgrid_); return; }

  // ================================================================================
  // Function needed to sort the SubGrids
  static bool ComparexMin(SubGrid const& sg1, SubGrid const& sg2)
  {
    if(sg1.xMin() == sg2.xMin()) throw runtime_error("There are SubGrids with the same lower bound.");
    return sg1.xMin() < sg2.xMin();
  }

  // ================================================================================
  // Create the global grid
  void grid::CreateGrid()
  {
    // Number of grids
    double ng = _GlobalGrid.size();

    // Check if there are extenal grids and if so disable the locking
    for(int ig=0; ig<ng; ig++) if(_GlobalGrid[ig].IsExt()) _ExtGrids = true;
    if(_ExtGrids && _Locked) {
      cout << endl;
      cout << "WARNING: External grids found." << endl;
      cout << "         ... unlocking SubGrids" << endl << endl;
      _Locked = false;
    }

    // Now oder the SubGrids in such a way that they start with that having the lowest value of xMin
    // (only if there is more than one grid).
    if(ng > 1) sort(_GlobalGrid.begin(), _GlobalGrid.end(), ComparexMin);

    // In case there grids have been locked ...
    if(_Locked) {

      // Find the point of the "(ig-1)"-th SubGrid such that "x[ig-1][ix] < xMin[ig] < x[ig-1][ix+1]",
      // and replace "xMin[ig]" with "x[ig-1][ix]"
      for(int ig=1; ig<ng; ig++) {

	int nxg     = _GlobalGrid[ig-1].nx();
	double xmin = _GlobalGrid[ig].xMin();

	// Parameters of the adjusted grid
	int nx_new;
	double xmin_new;
	int id_new = _GlobalGrid[ig].InterDegree();

	for(int ix=0; ix<nxg; ix++) {
	  double x = _GlobalGrid[ig-1].xg(ix);
	  if(x > xmin) {
	    xmin_new = x;
	    nx_new = nxg - ix;
	    break;
	  }
	}

	// Find the closest multiple of "nx - ix + 1" to "nx",
	// i.e. "DensityFactor", and replace "nx" accordingly.
	int DensityFactor = _GlobalGrid[ig].nx() / nx_new;
	nx_new *= DensityFactor;

	// Compute the new SubGrid and replace it in the global grid
	SubGrid sgrid(nx_new, xmin_new, id_new);
	_GlobalGrid[ig] = sgrid;
      }
    }

    // Compute the joint grid.
    // Parameters of the joint grid
    int nx_joint = -1;
    int id_joint = _GlobalGrid[0].InterDegree(); // Use the interpolation degree of the first grid
    vector<double> xg_joint_vect;

    double eps = 1e-12;
    for(int ig=0; ig<ng; ig++) {
      int nxg = _GlobalGrid[ig].nx();
      double xtrans;
      if(ig < ng-1) xtrans = _GlobalGrid[ig+1].xMin();
      else          xtrans = 1 + 2 * eps;
      for(int ix=0; ix<=nxg; ix++) {
	double x = _GlobalGrid[ig].xg(ix);

	if(xtrans - x < eps) break;
	nx_joint++;
	xg_joint_vect.push_back(x);
      }
    }

    // Copy the vector "xg_joint_vect" into the "xg_joint" and then
    // to initialize another SubGrid.
    double xg_joint[nx_joint];
    copy(xg_joint_vect.begin(), xg_joint_vect.end(), xg_joint);
    SubGrid grid_joint(nx_joint+1, xg_joint, id_joint);

    //// insert the joint SubGrid at the beginning of "_GlobalGrid"
    //_GlobalGrid.insert(_GlobalGrid.begin(), grid_joint);
    // Push the joint SubGrid at the end of "_GlobalGrid"
    _GlobalGrid.push_back(grid_joint);

    return;
  }

  // ================================================================================
  // Erase the global grid
  void grid::EraseGrid()
  {
    _GlobalGrid.erase(_GlobalGrid.begin(), _GlobalGrid.end());
    return;
  }

  int grid::nGrids() const
  {
    return _GlobalGrid.size() - 1;
  }

  // ================================================================================
  // Function to retrieve the ig-th subrid
  SubGrid const& grid::GetSubGrid(int ig) const
  {
    if(ig < 0 || ig >= _GlobalGrid.size()) throw runtime_error("Grid index out of range.");
    return _GlobalGrid[ig];
  }

  // ================================================================================
  // Check whether grids are equal
  bool grid::operator == (grid const& g)
  {
    vector<SubGrid> _GlobalGrid;
    if(_Locked != g._Locked)                       return false;
    if(_ExtGrids != g._ExtGrids)                   return false;
    if(_GlobalGrid.size() != g._GlobalGrid.size()) return false;
    for(int ig=0; ig<_GlobalGrid.size(); ig++) {
      if(!(_GlobalGrid[ig] == g._GlobalGrid[ig]))  return false;
    }
    return true;
  }

}
