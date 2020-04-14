//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include <vector>
#include <tuple>
#include <iostream>
#include <functional>

namespace apfel
{
  /**
   * @brief The template class QGrids is a mother class for the
   * interpolation in Q. This class also implements methods for the
   * subgrid interpolation relevant for example in a VFNS evolution.
   */
  template<class T>
  class QGrid
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    QGrid() = delete;

    /**
     * @brief The QGrid constructor.
     * @param nQ: the number of grid intervals in Q
     * @param QMin: the lower edge of the grid in Q
     * @param QMax: the upper edge of the grid in Q
     * @param InterDegree: the interpolation degree
     * @param Thresholds: the fixed point of the grid over which interpolation is forbidden
     * @param TabFunc: the function used to tabulate the grid in Q
     * @param InvTabFunc: the inverse function of TabFunc (an analytic expression is necessary)
     */
    QGrid(int                                  const& nQ,
          double                               const& QMin,
          double                               const& QMax,
          int                                  const& InterDegree,
          std::vector<double>                  const& Thresholds,
          std::function<double(double const&)> const& TabFunc,
          std::function<double(double const&)> const& InvTabFunc);

    /**
     * @brief The QGrid constructor.
     * @param nQ: the number of grid intervals in Q
     * @param QMin: the lower edge of the grid in Q
     * @param QMax: the upper edge of the grid in Q
     * @param InterDegree: the interpolation degree
     * @param Thresholds: the fixed point of the grid over which interpolation is forbidden
     * @param Lambda: the parameter of the function log(log(Q/Lambda)) used for the tabulation on the grid in Q
     */
    QGrid(int                 const& nQ,
          double              const& QMin,
          double              const& QMax,
          int                 const& InterDegree,
          std::vector<double> const& Thresholds,
          double              const& Lambda = 0.25);

    /**
     * @brief The QGrid constructor.
     * @param Qg: the user-defined interpolation grid
     * @param InterDegree: the interpolation degree
     * @note This constructor assumes that the function used to
     * compute the interpolating functions is linear. In addition, the
     * input vector 'Qg' is assumed to be ordered and
     * 'well-behaved'. To be used with care.
     */
    QGrid(std::vector<double> const& Qg,
          int                 const& InterDegree);
    ///@}

    /**
     * @brief Function that interpolates on the grid in Q.
     * @param Q: the value of the required interpolation
     * @return the interpolated value
     */
    T Evaluate(double const& Q) const;

    /**
     * @brief Function that derives on the grid in Q.
     * @param Q: the value of the required interpolation
     * @return the interpolated derivative
     */
    T Derive(double const& Q) const;

    /**
     * @brief Function that integrates on the grid in Q.
     * @param Qa: the lower integration bound
     * @param Qb: the upper integration bound
     * @return the interpolated integral
     */
    T Integrate(double const& Qa, double const& Qb) const;

    /**
     * @name Comparison operators
     * Collection of operators for comparing QGrid objects
     */
    ///@{
    bool operator == (QGrid const& sg) const;
    bool operator != (QGrid const& sg) const;
    ///@}

    /**
     * @name Getters
     */
    ///@{
    int                                         nQ()                 const { return _nQ; }          //!< return the number of Q interval
    int                                         InterDegree()        const { return _InterDegree; } //!< return the interpolation degree
    double                                      QMin()               const { return _QMin; }        //!< return the minimum node value
    double                                      QMax()               const { return _QMax; }        //!< return the maximum node value
    std::function<double(double const&)> const& TabFunc()            const { return _TabFunc; }     //!< return the tabulation function
    std::vector<double>                  const& GetThresholds()      const { return _Thresholds;}   //!< return the heavy quark thresholds
    std::vector<double>                  const& GetQGrid()           const { return _Qg; }          //!< return the grid in Q
    std::vector<double>                  const& GetFQGrid()          const { return _fQg; }         //!< return the grid in _TabFunc(Q)
    std::vector<int>                     const& GetThesholdIndices() const { return _nQg; }         //!< return the indices of the thresholds on the grid
    std::vector<T>                       const& GetQGridValues()     const { return _GridValues; }  //!< return the tabulated objects on the grid.
    ///@}

    /**
     * @brief Interpolation functions on QGrid.
     * @param tQ: interpolation control parameter
     * @param tau: the grid index
     * @param fq: the value of _TabFunc(Q) of the required interpolation
     * @return the interpolation weights
     */
    double Interpolant(int const& tQ, int const& tau, double const& fq) const;

    /**
     * @brief Derivative of the interpolation functions on QGrid.
     * @param tQ: interpolation control parameter
     * @param tau: the grid index
     * @param Q: the value of Q of the required interpolation
     * @return the derivarive of the interpolation weights
     */
    double DerInterpolant(int const& tQ, int const& tau, double const& Qa) const;

    /**
     * @brief Derivative of the interpolation functions on QGrid.
     * @param tQ: interpolation control parameter
     * @param tau: the grid index
     * @param Qa: the value of the lower integration bound
     * @param Qb: the value of the upper integration bound
     * @return the integral of interpolation weights
     */
    double IntInterpolant(int const& tQ, int const& tau, double const& Qa, double const& Qb) const;

    /**
     * @brief Computes the control parameter of the interpolant, the
     * lower and upper bounds over which the sum is limited.
     * @param Q: the value of the required interpolation
     * @return the lower and upper bounds of tau
     */
    std::tuple<int, int, int> SumBounds(double const& Q) const;

  protected:
    int                                  _nQ;           //!< Number intervals
    double                               _QMin;         //!< Minumim value of Q
    double                               _QMax;         //!< Maximum value of Q
    int                                  _InterDegree;  //!< Interpolation degree
    std::vector<double>                  _Thresholds;   //!< Thresholds
    std::function<double(double const&)> _TabFunc;      //!< Function whose constant step is used for the tabulation
    std::vector<double>                  _Qg;           //!< Grid in Q
    std::vector<double>                  _fQg;          //!< Grid in _TabFunc(Q)
    std::vector<int>                     _nQg;          //!< Indices of the nodes on which there is either a bound or a threshold
    std::vector<T>                       _GridValues;   //!< Vector of values to be interpolated on the grid

    template<class U>
    friend std::ostream& operator << (std::ostream& os, QGrid<U> const& dt);
  };

  /**
   * @brief Method that prints QGrid with cout <<.
   */
  template<class T>
  inline std::ostream& operator << (std::ostream& os, QGrid<T> const& Qg)
  {
    os << "QGrid: " << &Qg << "\n";
    os << "nQ                = " << Qg._nQ << "\n";
    os << "QMin              = " << Qg._QMin << "\n";
    os << "QMax              = " << Qg._QMax << "\n";
    os << "InterDegree       = " << Qg._InterDegree << "\n";
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
