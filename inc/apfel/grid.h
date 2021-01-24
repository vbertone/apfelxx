//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/subgrid.h"

#include <memory>

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
     * @brief The Grid constructor.
     * @param grs: vector of subgrids
     */
    Grid(std::vector<SubGrid> const& grs);

    /**
     * @name Getters
     */
    ///@{
    /**
     * @return The number of subgrids
     */
    int nGrids() const { return _GlobalGrid.size(); }

    /**
     * @return The map of indices from the joint grid to the subgrids
     */
    std::vector<std::pair<int, int>> const& SubToJointMap() const { return _SubToJointMap; }

    /**
     * @return The map of indices from the subgrids to the joint grid
     */
    std::vector<std::vector<int>> const& JointToSubMap() const { return _JointToSubMap; }

    /**
     * @return The vector of transition indices on the joint grid
     */
    std::vector<int> const& TransitionPoints() const { return _TransPoints; }

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
    SubGrid const& GetJointGrid() const { return *_JointGrid; }
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
     * @brief Fill in the joint grid object with the appropriate grid
     * nodes.
     * @return the joint grid
     */
    SubGrid CreateJointGrid();

  private:
    std::vector<std::pair<int, int>> _SubToJointMap; //!< Vector of pairs corresponding to grid- and node-indices on the subgrids.
    std::vector<std::vector<int>>    _JointToSubMap; //!< Vector of indices from the subgrids to the joint grid
    std::vector<int>                 _TransPoints;   //!< Vector of indices corresponding to the transition from one subgrid to the other
    std::vector<SubGrid>             _GlobalGrid;    //!< Vector with sub-grids.
    std::unique_ptr<SubGrid>         _JointGrid;     //!< Container for the joint grid.

    friend std::ostream& operator << (std::ostream& os, Grid const& gr);
  };

  /**
   * @brief Overload the << operator to print the parameters of the
   * grid.
   */
  std::ostream& operator << (std::ostream& os, Grid const& gr);
}
