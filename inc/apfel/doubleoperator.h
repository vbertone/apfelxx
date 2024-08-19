//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/doubleexpression.h"
#include "apfel/grid.h"
#include "apfel/matrix.h"
#include "apfel/operator.h"
#include "apfel/doubledistribution.h"

namespace apfel
{
  /**
   * @brief The DoubleOperator class defines is essentially the
   * convolution on a pair of grids bewteen an DoubleExpression object
   * and the interpolant functions.
   */
  class DoubleOperator
  {
  public:
    DoubleOperator() = delete;
    DoubleOperator(DoubleOperator const&) = default;

    /**
     * @brief The DoubleOperator constructor.
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param dexpr: the double expression to be transformed
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     */
    DoubleOperator(Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr, double const& eps = 1e-5);

    /**
     * @brief The DoubleOperator constructor.
     * @param O1: the first single operator
     * @param O2: the second single operator
     * @param dexpr: the double expression to be transformed (Null by default)
     */
    DoubleOperator(Operator const& O1, Operator const& O2, DoubleExpression const& dexpr = DoubleExpression{});

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~DoubleOperator() {}

    /**
     * @name Binary operators
     */
    ///@{
    DoubleDistribution operator *= (DoubleDistribution const& d) const;                  //!< this *= Distribution
    DoubleOperator& operator *= (double const& s);                                       //!< this *= Scalar
    DoubleOperator& operator /= (double const& s);                                       //!< this /= Scalar
    DoubleOperator& operator += (DoubleOperator const& o);                               //!< this += Operator
    DoubleOperator& operator -= (DoubleOperator const& o);                               //!< this -= Operator
    DoubleOperator& operator *= (DoubleOperator const& o);                               //!< this *= Operator
    DoubleOperator& operator  = (DoubleOperator const& o);                               //!< this  = Operator
    DoubleOperator& operator *= (std::function<double(double const&, double const&)> f); //!< This *= Function
    ///@}

    /**
     * @brief Function that returns the first Grid object associated
     * to the double operator.
     */
    Grid const& GetFirstGrid() const { return _grid1; }

    /**
     * @brief Function that returns the second Grid object associated
     * to the double operator.
     */
    Grid const& GetSecondGrid() const { return _grid2; }

    /**
     * @brief Function that returns the DoubleExpression object
     * associated to the operator.
     */
    DoubleExpression const& GetDoubleExpression() const { return _dexpr; }

    /**
     * @brief Function that returns the integration accuracy
     */
    double const& GetIntegrationAccuracy() const { return _eps; }

    /**
     * @brief Function that returns the double operator.
     */
    std::vector<std::vector<matrix<matrix<double>>>> GetDoubleOperator() const { return _dOperator; }

    /**
     * @brief Print the Operator object
     */
    void Print() const { std::cout << *this << std::endl; }

  protected:
    Grid                                      const& _grid1;      //!< First grid on which to compute the operator
    Grid                                      const& _grid2;      //!< Second grid on which to compute the operator
    DoubleExpression                          const& _dexpr;      //!< Expression to be used
    double                                    const  _eps;        //!< Integration accuracy
    std::vector<std::vector<matrix<matrix<double>>>> _dOperator;  //!< DoubleOperator values

    friend std::ostream& operator << (std::ostream& os, DoubleOperator const& dop);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution operator * (DoubleOperator lhs, DoubleDistribution const& rhs);                        //!< Operator*Distribution
  DoubleOperator     operator * (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< Operator*Operator
  DoubleOperator     operator * (double const& s, DoubleOperator rhs);                                      //!< Scalar*Operator
  DoubleOperator     operator * (DoubleOperator lhs, double const& s);                                      //!< Operator*Scalar
  DoubleOperator     operator * (std::function<double(double const&, double const)> f, DoubleOperator rhs); //!< function*Operator
  DoubleOperator     operator * (DoubleOperator lhs, std::function<double(double const&, double const)> f); //!< Operator*function
  DoubleOperator     operator / (DoubleOperator lhs, double const& s);                                      //!< Operator/Scalar
  DoubleOperator     operator + (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< Operator+Operator
  DoubleOperator     operator - (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< Operator-Operator
  ///@}

  /**
   * @brief Method which prints DoubleOperator with std::cout <<.
   */
  std::ostream& operator << (std::ostream& os, DoubleOperator const& dop);
}

