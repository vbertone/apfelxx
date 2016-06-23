/*
  grid.hh:

  Author: Valerio Bertone
*/

#pragma once

#include <vector>

#include "subgrid.hh"

using namespace std;

namespace apfel {

  /**
   * class grid:
   * class for the global x-grid object.
   * This class defines the "grid" object which is essentially a collection of
   * "subgrid" objects plus other global parameters. This class also includes
   * all the relevant methods for the manipulation of the subgrids.
   */
  class grid {
  private:
    // Attributes
    vector <subgrid> _GlobalGrid;
    bool _Locked;
    bool _ExtGrids;

  public:
    // Constructors
    grid();

    // Destructor
    ~grid();

    // Setters
    void LockGrids() { _Locked = true; }
    void AddGrid(subgrid const& sgrid_);
    void CreateGrid();
    void EraseGrid();

    // Getters
    int     nGrids()        const { return _GlobalGrid.size() - 1; }
    bool    Locked()        const { return _Locked; }
    bool    ExtGrids()      const { return _ExtGrids; }
    subgrid SubGrid(int ig) const;

    // Check whether grids are equal
    bool operator == (grid const& g);
  };

}
