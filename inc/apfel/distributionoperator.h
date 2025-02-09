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
   * @brief The DistributionOperator class defines an object that
   * behaves as a distribution along the first variable and as an
   * operator along the second.
   */
  class DistributionOperator
  {
  public:
    DistributionOperator() = delete;
    DistributionOperator(DistributionOperator const&) = default;

    /**
     * @brief The DistributionOperator constructor.
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     */
    DistributionOperator(Grid const& gr1, Grid const& gr2);

    /**
     * @brief The DistributionOperator constructor.
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param DO: the DistributionOperator container
     */
    DistributionOperator(Grid const& gr1, Grid const& gr2, std::vector<std::vector<std::vector<matrix<double>>>> const& DO);

    /**
     * @brief The DistributionOperator constructor.
     * @param d1: the distribution on the first variable
     * @param O2: the operator on the second variable
     */
    DistributionOperator(Distribution const& d1, Operator const& O2);

    /**
     * @brief The DistributionOperator constructor.
     * @param DObj: DoubleObject of distributions and operators
     */
    DistributionOperator(DoubleObject<Distribution, Operator> const& DObj);

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~DistributionOperator() {}

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
    DistributionOperator& operator *= (Operator const& o);                                            //!< this *= Operator
    DistributionOperator& operator *= (double const& s);                                              //!< this *= Scalar
    DistributionOperator& operator /= (double const& s);                                              //!< this /= Scalar
    DistributionOperator& operator += (DistributionOperator const& o);                                //!< this += DistributionOperator
    DistributionOperator& operator -= (DistributionOperator const& o);                                //!< this -= DistributionOperator
    DistributionOperator& operator  = (DistributionOperator const& o);                                //!< this  = DistributionOperator
    DistributionOperator& operator *= (std::function<double(double const&, double const&)> const& f); //!< this *= Function of both variables
    DistributionOperator& operator *= (std::function<double(double const&)> const& f);                //!< this *= Function of one variable used for both variables
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
     * @brief Function that returns the DistributionOperator container.
     */
    std::vector<std::vector<std::vector<matrix<double>>>> GetDistributionOperator() const { return _dOperator; }

    /**
     * @brief Function that prints the DistributionOperator object
     */
    void Print() const { std::cout << *this << std::endl; }

  private:
    Grid                                           const& _grid1;      //!< First grid
    Grid                                           const& _grid2;      //!< Second grid
    LagrangeInterpolator                           const  _li1;        //!< The Lagrange interpolator on the first grid
    std::vector<std::vector<std::vector<matrix<double>>>> _dOperator;  //!< DistributionOperator container

    friend std::ostream& operator << (std::ostream& os, DistributionOperator const& dop);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution   operator * (DistributionOperator const& lhs, Distribution const& rhs);                               //!< DistributionOperator*Distribution
  DistributionOperator operator * (DistributionOperator lhs, Operator const& rhs);                                          //!< DistributionOperator*Operator
  DistributionOperator operator * (double const& s, DistributionOperator rhs);                                              //!< Scalar*DistributionOperator
  DistributionOperator operator * (DistributionOperator lhs, double const& s);                                              //!< DistributionOperator*Scalar
  DistributionOperator operator / (DistributionOperator lhs, double const& s);                                              //!< DistributionOperator/Scalar
  DistributionOperator operator + (DistributionOperator lhs, DistributionOperator const& rhs);                              //!< DistributionOperator+Operator
  DistributionOperator operator - (DistributionOperator lhs, DistributionOperator const& rhs);                              //!< DistributionOperator-DistributionOperator
  DistributionOperator operator * (std::function<double(double const&, double const&)> const& f, DistributionOperator rhs); //!< Function*DistributionOperator
  DistributionOperator operator * (DistributionOperator lhs, std::function<double(double const&, double const&)> const& f); //!< DistributionOperator*Function
  DistributionOperator operator * (std::function<double(double const&)> const& f, DistributionOperator rhs);                //!< Function*DistributionOperator
  DistributionOperator operator * (DistributionOperator lhs, std::function<double(double const&)> const& f);                //!< DistributionOperator*Function
  ///@}

  /**
   * @brief Method which prints DistributionOperator with cout <<.
   */
  std::ostream& operator << (std::ostream& os, DistributionOperator const& in);
}

