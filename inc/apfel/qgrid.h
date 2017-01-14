//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>
using std::vector;
#include <tuple>
using std::tuple;
using std::get;
#include <iostream>

namespace apfel
{
  /**
   * @brief Class for the Q-space interpolation QGrids.
   *
   * Subgrids are the building blocks of the interpolation procedure.
   * This class defines the "QGrid" object that includes, apart from the
   * grid itself, also the relevant parameters.
   *
   * This template class accepts fixed types:
   * - double
   *
   */
  template<class T>
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
     * @brief Interpolate on the grid
     * @param Q the value of the required interpolation
     * @return the interpolated value.
     */
    double Evaluate(double const& Q) const;

    /**
     * @brief Check whether QGrids are equal
     * @param sg the QGrid to be compared
     * @return true/false
     */
    bool operator == (QGrid const& sg) const;
    bool operator != (QGrid const& sg) const;

    // Getters
    int                   nQ()            const { return _nQ; }          //!< return the number of Q interval
    int                   InterDegree()   const { return _InterDegree; } //!< return the interpolation degree
    double                QMin()          const { return _QMin; }        //!< return the minimum node value
    double                QMax()          const { return _QMax; }        //!< return the maximum node value
    double                Lambda()        const { return _Lambda; }      //!< return the the value of _Lambda
    vector<double> const& GetThresholds() const { return _Thresholds;}   //!< return the heavy quark thresholds
    vector<double> const& GetQGrid()      const { return _Qg; }          //!< return the grid in Q
    vector<double> const& GetLogQGrid()   const { return _llQ2g; }       //!< return the grid in ln(ln(Q2/Lambda))

  protected:
    /**
     * @brief Interpolation functions on QGrid
     * @param tQ interpolation control parameter
     * @param tau the grid index
     * @param lnln2ql the value of ln(ln(Q^2/Lambda^2)) of the required interpolation
     * @return the interpolation weights.
     */
    double Interpolant(int const& tQ, int const& tau, double const& lnln2ql) const;

    /**
     * @brief Computes the control parameter of the interpolant, the lower and upper bounds over which the sum is limited
     * @param Q the value of the required interpolation
     * @return the lower and upper bounds of tau.
     */
    tuple<int,int,int> SumBounds(double const& Q) const;

    int            _nQ;           //!< Number intervals
    int            _InterDegree;  //!< Interpolation degree
    double         _QMin;         //!< Minumim value of Q
    double         _QMax;         //!< Maximum value of Q
    double         _Lambda;       //!< Value of the pole of the distribution of the nodes
    vector<double> _Thresholds;   //!< Heavy quark thresholds
    vector<double> _Qg;           //!< Actual grid in Q
    vector<double> _llQ2g;        //!< Actual grid in ln(ln(Q^2/Lambda^2))
    vector<int>    _nQg;          //!< Indices of the nodes on which there is either a bound or a threshold

    vector<T> _GridValues;   //!< Vector of values to be interpolated on the grid

    template<class U>
    friend std::ostream& operator<<(std::ostream& os, QGrid<U> const& dt);
  };

  /**
   * @brief Method which prints QGrid with cout <<.
   */
  template<class T>
  std::ostream& operator<<(std::ostream& os, QGrid<T> const& Qg)
  {
    os << "QGrid: " << &Qg << "\n";
    os << "nQ                = " << Qg._nQ << "\n";
    os << "QMin              = " << Qg._QMin << "\n";
    os << "QMax              = " << Qg._QMax << "\n";
    os << "InterDegree       = " << Qg._InterDegree << "\n";
    os << "Lambda            = " << Qg._Lambda << "\n";
    os << "Thresholds        = ";
    for (const auto &v: Qg._Thresholds) os << v << " ";
    os << "\n";
    os << "Threshold indices = ";
    for (const auto &v: Qg._nQg) os << v << " ";
    os << "\n";
    os << "Qg                = ";
    for (const auto &v: Qg._Qg) os << v << " ";
    os << "\n";
    return os;
  }
}
