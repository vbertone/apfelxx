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

namespace apfel
{
  /**
   * @brief The DistributionOperator class defines an object that
   * behaves as a distribution for the first variable and as an
   * operator for the second.
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
     * @param d1: the distribution on first
     * @param O2: the operator on the second variable
     */
    DistributionOperator(Distribution const& d1, Operator const& O2);

    /**
     * @brief The DistributionOperator constructor.
     * @param DObj: double object of operators
     */
    DistributionOperator(DoubleObject<Distribution, Operator> const& DObj);

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~DistributionOperator() {}

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
    DistributionOperator& operator *= (std::function<double(double const&, double const&)> const& f); //!< this *= Function of both the integration variables
    DistributionOperator& operator *= (std::function<double(double const&)> const& f);                //!< this *= Function of one variable used for both the integration variables
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
     * @brief Function that returns the double operator.
     */
    std::vector<std::vector<std::vector<matrix<double>>>> GetDistributionOperator() const { return _dOperator; }

    /**
     * @brief Print the DistributionOperator object
     */
    void Print() const { std::cout << *this << std::endl; }

  private:
    Grid                                           const& _grid1;      //!< First grid on which to compute the operator
    Grid                                           const& _grid2;      //!< Second grid on which to compute the operator
    std::vector<std::vector<std::vector<matrix<double>>>> _dOperator;  //!< DistributionOperator values

    friend std::ostream& operator << (std::ostream& os, DistributionOperator const& dop);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution   operator * (DistributionOperator lhs, Distribution const& rhs);                                      //!< DistributionOperator*Distribution
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

