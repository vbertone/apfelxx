//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>
using std::vector;

namespace apfel
{
  /**
   * @brief Class for the Q-space interpolation QGrids.
   *
   * Subgrids are the building blocks of the interpolation procedure.
   * This class defines the "QGrid" object that includes, apart from the
   * grid itself, also the relevant parameters.
   *
   */
  class QGrid
  {
  public:

    QGrid() = delete;

    /**
     * @brief Standard internal grid constructor.
     * @param nQ number of grid intervals in Q.
     * @param QMin lower edge Q of the grid.
     * @param QMax upper edge Q of the grid.
     * @param InterDegree interpolation degree.
     */
    QGrid(int const& nQ, double const& QMin, double const& QMax, int const& InterDegree, vector<double> const& Thresholds, double const& Lambda = 0.25);

    /**
     * @brief Check whether QGrids are equal
     * @param sg the QGrid to be compared
     * @return true/false
     */
    bool operator == (QGrid const& sg) const;
    bool operator != (QGrid const& sg) const;

    // Getters
    int                   nQ()          const { return _nQ; }          //!< return the number of Q interval
    int                   InterDegree() const { return _InterDegree; } //!< return the interpolation degree
    double                QMin()        const { return _QMin; }        //!< return the minimum node value
    double                QMax()        const { return _QMax; }        //!< return the maximum node value
    vector<double> const& GetQGrid()    const { return _Qg; }          //!< return the grid in Q
    vector<double> const& GetLogQGrid() const { return _llQ2g; }       //!< return the grid in ln(ln(Q2/Lambda))

  private:
    int            _nQ;           //!< Number intervals
    int            _InterDegree;  //!< Interpolation degree
    double         _QMin;         //!< Minumim value of Q
    double         _QMax;         //!< Maximum value of Q
    vector<double> _Thresholds;   //!< Heavy quark thresholds
    double         _Lambda;       //!< Value of the pole of the distribution of the nodes
    vector<double> _Qg;           //!< Actual grid in Q
    vector<double> _llQ2g;        //!< Actual grid in ln(ln(Q^2/Lambda^2))
    vector<int>    _nQg;          //!< Indices of the nodes on which there is either a bound or a threshold

    friend std::ostream& operator<<(std::ostream& os, QGrid const& dt);
  };

  /**
   * @brief Method which prints QGrid with cout <<.
   */
  std::ostream& operator<<(std::ostream& os, QGrid const& sg);
}
