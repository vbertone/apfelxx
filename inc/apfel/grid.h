/*
  grid.hh:

  Author: Valerio Bertone
*/

#pragma once

#include <vector>
#include <apfel/subgrid.h>
using std::vector;

namespace apfel {

  /**
   * class grid:
   * class for the global x-grid object.
   * This class defines the "grid" object which is essentially a collection of
   * "SubGrid" objects plus other global parameters. This class also includes
   * all the relevant methods for the manipulation of the SubGrids.
   */
  class grid
  {
  private:
    // Attributes
    vector<SubGrid> _GlobalGrid;
    bool _Locked;
    bool _ExtGrids;

  public:
    // Constructors
    grid();

    // Destructor
    ~grid();

    // Setters
    void LockGrids() { _Locked = true; }
    void AddGrid(SubGrid const& sgrid_);
    void CreateGrid();
    void EraseGrid();

    // Getters
    int     nGrids() const;
    bool    Locked()        const { return _Locked; }
    bool    ExtGrids()      const { return _ExtGrids; }
    SubGrid const& GetSubGrid(int ig) const;

    // Check whether grids are equal
    bool operator == (grid const& g);
  };

}
