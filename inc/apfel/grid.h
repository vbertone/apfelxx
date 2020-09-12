//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/subgrid.h"

namespace apfel
{
  /**
   * @brief The Grid class defines ab object that is essentially a
   * collection of "SubGrid" objects plus other global
   * parameters. This class also includes all the relevant methods for
   * the manipulation of the SubGrids.
   */
  class Grid
  {
  public:
    /**
     * @brief Default constructor.
     */
    Grid();

    /**
     * @brief The Grid constructor.
     * @param grs: vector of subgrids
     * @param lockgrids: flag to enable/disable the locking of the subgrids (default: true)
     */
    Grid(std::vector<SubGrid> const& grs, bool const& lockgrids = true);

    /**
     * @name Getters
     */
    ///@{
    /**
     * @return The number of subgrids
     */
    int nGrids() const { return _GlobalGrid.size(); }

    /**
     * @return The locking flag
     */
    bool Locked() const { return _Locked; }

    /**
     * @return The external-grid flag
     */
    bool ExtGrids() const { return _ExtGrids; }

    /**
     * @return The map of indices from the subgrids to the joint grid
     */
    std::vector<std::vector<int>> const& JointToSubMap() const { return _JointToSubMap; }

    /**
     * @return The vector of subgrids
     */
    std::vector<SubGrid> const& GetSubGrids() const { return _GlobalGrid; }

    /**
     * @return The ig-th SubGrid
     */
    SubGrid const& GetSubGrid(int ig) const { return _GlobalGrid[ig]; }

    /**
     * @return The joint SubGrid
     */
    SubGrid const& GetJointGrid() const { return _JointGrid; }
    ///@}

    /**
     * @name Comparison operators
     * Collection of operators for comparing grid objects
     */
    ///@{
    bool operator == (Grid const& g) const;
    bool operator != (Grid const& g) const;
    ///@}

  private:
    /**
     * @brief Takes the input SubGrids, apply the locking if needed
     * and fill the joint grid object with the appropriate grid nodes.
     * @return the joint grid
     */
    SubGrid CreateJointGrid();

  private:
    bool                          _Locked;        //!< Flag for locking the grids.
    bool                          _ExtGrids;      //!< Contains external sub-grids.
    std::vector<std::vector<int>> _JointToSubMap; //!< Map of indices from the subgrids to the joint grid (empty if the subgrids are not locked)
    std::vector<SubGrid>          _GlobalGrid;    //!< Vector with sub-grids.
    SubGrid                       _JointGrid;     //!< Container for the joint grid.

    friend std::ostream& operator << (std::ostream& os, Grid const& gr);
  };

  /**
   * @brief Overload the << operator to print the parameters of the
   * grid.
   */
  std::ostream& operator << (std::ostream& os, Grid const& gr);
}
