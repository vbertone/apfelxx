//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/matrix.h"
#include "apfel/operator.h"
#include "apfel/distribution.h"
#include "apfel/doubleobject.h"
#include "apfel/doubledistribution.h"
#include "apfel/lagrangeinterpolator.h"

namespace apfel
{
  /**
   * @brief The OperatorDistribution class defines an object that
   * behaves as an operator along the first variable and as a
   * distribution along the second.
   */
  class OperatorDistribution
  {
  public:
    OperatorDistribution() = delete;
    OperatorDistribution(OperatorDistribution const&) = default;

    /**
     * @brief The OperatorDistribution constructor.
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     */
    OperatorDistribution(Grid const& gr1, Grid const& gr2);

    /**
     * @brief The OperatorDistribution constructor.
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param OD: the OperatorDistribution container
     */
    OperatorDistribution(Grid const& gr1, Grid const& gr2, std::vector<std::vector<matrix<std::vector<double>>>> const& OD);

    /**
     * @brief The OperatorDistribution constructor.
     * @param O1: the operator on the first variable
     * @param d2: the distribution on the second variable
     */
    OperatorDistribution(Operator const& O1, Distribution const& d2);

    /**
     * @brief The OperatorDistribution constructor.
     * @param DObj: DoubleObject of operators and distributions
     */
    OperatorDistribution(DoubleObject<Operator, Distribution> const& DObj);

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~OperatorDistribution() {}

    /**
     * @name Evaluate functions
     * List of functions that perform the interpolation on the x-space
     * grids.
     */
    ///@{
    /**
     * @brief Function that evaluates the Distribution part in a given
     * point through interpolation.
     * @param x: the value in x where the distribution part is to be interpolated
     * @return the interpolated result
     */
    Operator Evaluate(double const& x) const;

    /**
     * @brief Function that evaluates the derivative of the Distribution
     * part in a given point through interpolation.
     * @param x: the value in x where the distribution part is to be interpolated
     * @return the derivative of the interpolated function
     */
    Operator Derive(double const& x) const;

    /**
     * @brief Function that evaluates the integral of the Distribution
     * part in the interval[a, b].
     * @param a: the lower integration bound for the first variable
     * @param b: the upper integration bound for the first variable
     * @return the integral of the interpolated function
     */
    Operator Integrate(double const& a, double const& b) const;
    ///@}

    /**
     * @name Binary operators
     */
    ///@{
    DoubleDistribution operator    *= (Distribution const& d) const;                                  //!< this *= Distribution
    OperatorDistribution& operator *= (Operator const& o);                                            //!< this *= Operator
    OperatorDistribution& operator *= (double const& s);                                              //!< this *= Scalar
    OperatorDistribution& operator /= (double const& s);                                              //!< this /= Scalar
    OperatorDistribution& operator += (OperatorDistribution const& o);                                //!< this += OperatorDistribution
    OperatorDistribution& operator -= (OperatorDistribution const& o);                                //!< this -= OperatorDistribution
    OperatorDistribution& operator  = (OperatorDistribution const& o);                                //!< this  = OperatorDistribution
    OperatorDistribution& operator *= (std::function<double(double const&, double const&)> const& f); //!< this *= Function of both variables
    OperatorDistribution& operator *= (std::function<double(double const&)> const& f);                //!< this *= Function of one variable used for both variables
    ///@}

    /**
     * @brief Function that returns the first Grid object.
     */
    Grid const& GetFirstGrid() const { return _grid1; }

    /**
     * @brief Function that returns the second Grid object.
     */
    Grid const& GetSecondGrid() const { return _grid2; }

    /**
     * @brief Function that returns the OperatorDistribution container.
     */
    std::vector<std::vector<matrix<std::vector<double>>>> GetOperatorDistribution() const { return _dOperator; }

    /**
     * @brief Function that prints the OperatorDistribution object
     */
    void Print() const { std::cout << *this << std::endl; }

  private:
    Grid                                           const& _grid1;      //!< First grid
    Grid                                           const& _grid2;      //!< Second grid
    LagrangeInterpolator                           const  _li2;        //!< The Lagrange interpolator on the second grid
    std::vector<std::vector<matrix<std::vector<double>>>> _dOperator;  //!< OperatorDistribution container

    friend std::ostream& operator << (std::ostream& os, OperatorDistribution const& opd);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution   operator * (OperatorDistribution const& lhs, Distribution const& rhs);                               //!< OperatorDistribution*Distribution
  OperatorDistribution operator * (OperatorDistribution lhs, Operator const& rhs);                                          //!< OperatorDistribution*Operator
  OperatorDistribution operator * (double const& s, OperatorDistribution rhs);                                              //!< Scalar*OperatorDistribution
  OperatorDistribution operator * (OperatorDistribution lhs, double const& s);                                              //!< OperatorDistribution*Scalar
  OperatorDistribution operator / (OperatorDistribution lhs, double const& s);                                              //!< OperatorDistribution/Scalar
  OperatorDistribution operator + (OperatorDistribution lhs, OperatorDistribution const& rhs);                              //!< OperatorDistribution+Operator
  OperatorDistribution operator - (OperatorDistribution lhs, OperatorDistribution const& rhs);                              //!< OperatorDistribution-OperatorDistribution
  OperatorDistribution operator * (std::function<double(double const&, double const&)> const& f, OperatorDistribution rhs); //!< Function*OperatorDistribution
  OperatorDistribution operator * (OperatorDistribution lhs, std::function<double(double const&, double const&)> const& f); //!< OperatorDistribution*Function
  OperatorDistribution operator * (std::function<double(double const&)> const& f, OperatorDistribution rhs);                //!< Function*OperatorDistribution
  OperatorDistribution operator * (OperatorDistribution lhs, std::function<double(double const&)> const& f);                //!< OperatorDistribution*Function
  ///@}

  /**
   * @brief Method which prints OperatorDistribution with cout <<.
   */
  std::ostream& operator << (std::ostream& os, OperatorDistribution const& in);
}

