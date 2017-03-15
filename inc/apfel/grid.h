//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/subgrid.h>

#include <vector>
#include <memory>

using std::vector;
using std::unique_ptr;

namespace apfel
{
  /**
   * @brief The global x-grid object.
   *
   * This class defines the "Grid" object which is essentially a collection of
   * "SubGrid" objects plus other global parameters. This class also includes
   * all the relevant methods for the manipulation of the SubGrids.
   */
  class Grid
  {
  public:

    Grid() = delete;

    /**
     * @brief The default constructor, takes a vector of SubGrid
     * @param grs vector of subgrids
     * @param lockgrids flag to enable-disable the locking.
     */
    Grid(vector<SubGrid> const& grs, bool const& lockgrids = true);

    // Getters
    int            nGrids()           const { return _GlobalGrid.size(); } //!< return the number of subgrids
    bool           Locked()           const { return _Locked; }            //!< return the locked flag
    bool           ExtGrids()         const { return _ExtGrids; }          //!< return the external grid flag
    SubGrid const& GetSubGrid(int ig) const { return _GlobalGrid[ig]; }    //!< return the SubGrid item ig
    SubGrid const& GetJointGrid()     const { return *_JointGrid; }        //!< return the joint grid

    /**
     * @brief Check whether Grids are equal
     * @param sg the SubGrid to be compared
     * @return true/false
     */
    bool operator == (Grid const& g) const;
    bool operator != (Grid const& g) const;

  private:
    /**
     * @brief Takes the input SubGrids, apply the locking and
     * fill the \c _JointGrid object with the appropriate grid nodes.
     */
    void CreateJointGrid();

    bool                _Locked;     //!< Flag for locking the grids.
    bool                _ExtGrids;   //!< Contains external sub-grids.
    vector<SubGrid>     _GlobalGrid; //!< Vector with sub-grids.
    unique_ptr<SubGrid> _JointGrid;  //!< Container for the joint grid.

    friend std::ostream& operator<<(std::ostream& os, Grid const& gr);
  };

  /**
   * @brief Method which prints Grid with cout <<.
   */
  std::ostream& operator<<(std::ostream& os, Grid const& gr);

}
