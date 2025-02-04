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
#include "apfel/distributionoperator.h"
#include "apfel/operatordistribution.h"
#include "apfel/config.h"

// Include yaml-cpp header only if it has been found at configuration
// time.
#if WITH_YAML_CPP == 1
#include <yaml-cpp/yaml.h>
#endif

namespace apfel
{
  /**
   * @brief This class defines the basic object DoubleOperator which
   * essentially contains the convolution bewteen a "DoubleExpression"
   * object and two sets of interpolanting functions defined on two
   * independent grids. This is the two-dimensional generalisation of
   * the object "Operator".
   * @note For the moment, only the inclusive massless case can be
   * handled within this class.
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
     * @param dexpr: the double expression to be convoluted
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-3</SUP>)
     */
    DoubleOperator(Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr, double const& eps = 1e-3);

    /**
     * @brief The DoubleOperator constructor.
     * @param O1: the first single operator
     * @param O2: the second single operator
     */
    DoubleOperator(Operator const& O1, Operator const& O2);

    /**
     * @brief The DoubleOperator constructor.
     * @param DObj: double object of operators
     */
    DoubleOperator(DoubleObject<Operator> const& DObj);

#if WITH_YAML_CPP == 1
    /**
     * @brief The DoubleOperator constructor.
     * @param Node: YAML Node where the DoubleOperator object is strored
     * @param gr1: the Grid object for the first variable
     * @param gr2: the Grid object for the second variable
     * @param dexpr: the double expression to be convoluted needed to compare names
     */
    DoubleOperator(YAML::Node const& Node, Grid const& gr1, Grid const& gr2, DoubleExpression const& dexpr);

    /**
     * @brief Function the encapsulate a DoubleOperator object in a
     * YAML::Emitter objects such that it can be written on file.
     */
    std::string EmitDoubleOperator() const;
#endif

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~DoubleOperator() {}

    /**
     * @name Binary operators
     */
    ///@{
    DoubleDistribution operator *= (DoubleDistribution const& d) const;                     //!< this *= DoubleDistribution
    DistributionOperator MultiplyFirstBy(Distribution const& d) const;                      //!< this times a distribution convoluted with the first variable
    OperatorDistribution MultiplySecondBy(Distribution const& d) const;                     //!< this times a distribution convoluted with the second variable
    DoubleOperator&    operator *= (double const& s);                                       //!< this *= Scalar
    DoubleOperator&    operator /= (double const& s);                                       //!< this /= Scalar
    DoubleOperator&    operator += (DoubleOperator const& o);                               //!< this += DoubleOperator
    DoubleOperator&    operator -= (DoubleOperator const& o);                               //!< this -= DoubleOperator
    DoubleOperator&    operator *= (DoubleOperator const& o);                               //!< this *= DoubleOperator
    DoubleOperator&    operator  = (DoubleOperator const& o);                               //!< this  = DoubleOperator
    DoubleOperator&    operator *= (std::function<double(double const&, double const&)> f); //!< this *= 2D-function
    DoubleOperator&    operator *= (std::function<double(double const&)> f);                //!< this *= 1D-Function
    ///@}

    /**
     * @brief Function that returns the first Grid object associated
     * with the double operator.
     */
    Grid const& GetFirstGrid() const { return _grid1; }

    /**
     * @brief Function that returns the second Grid object associated
     * with the double operator.
     */
    Grid const& GetSecondGrid() const { return _grid2; }

    /**
     * @brief Function that returns the integration accuracy.
     */
    double const& GetIntegrationAccuracy() const { return _eps; }

    /**
     * @brief Function that returns the name of the DoubleExpression
     * object associated with the operator.
     */
    std::string const& GetDoubleExpressionName() const { return _dexprName; }

    /**
     * @brief Function that returns the DoubleOperator container.
     */
    std::vector<std::vector<matrix<matrix<double>>>> GetDoubleOperator() const { return _dOperator; }

    /**
     * @brief Function that prints the DoubleOperator object.
     */
    void Print() const { std::cout << *this << std::endl; }

  protected:
    Grid                                      const& _grid1;      //!< First grid on which to compute the operator
    Grid                                      const& _grid2;      //!< Second grid on which to compute the operator
    double                                    const  _eps;        //!< Integration accuracy
    std::string                               const  _dexprName;  //!< Name of the double expression
    std::vector<std::vector<matrix<matrix<double>>>> _dOperator;  //!< DoubleOperator container

    friend std::ostream& operator << (std::ostream& os, DoubleOperator const& dop);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  DoubleDistribution operator * (DoubleOperator const& lhs, DoubleDistribution const& rhs);                 //!< DoubleOperator*DoubleDistribution
  DistributionOperator operator * (Distribution const& lhs, DoubleOperator const& rhs);                     //!< Distribution*DoubleOperator = DistributionOperator
  OperatorDistribution operator * (DoubleOperator const& lhs, Distribution const& rhs);                     //!< DoubleOperator*Distribution = OperatorDistribution
  DoubleOperator     operator * (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< DoubleOperator*DoubleOperator
  DoubleOperator     operator * (double const& s, DoubleOperator rhs);                                      //!< Scalar*DoubleOperator
  DoubleOperator     operator * (DoubleOperator lhs, double const& s);                                      //!< DoubleOperator*Scalar
  DoubleOperator     operator * (std::function<double(double const&, double const)> f, DoubleOperator rhs); //!< 2D-function*DoubleOperator
  DoubleOperator     operator * (DoubleOperator lhs, std::function<double(double const&, double const)> f); //!< DoubleOperator*2D-function
  DoubleOperator     operator * (std::function<double(double const&)> f, DoubleOperator rhs);               //!< 1D-function*DoubleOperator
  DoubleOperator     operator * (DoubleOperator lhs, std::function<double(double const&)> f);               //!< DoubleOperator*1D-function
  DoubleOperator     operator / (DoubleOperator lhs, double const& s);                                      //!< DoubleOperator/Scalar
  DoubleOperator     operator + (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< DoubleOperator+DoubleOperator
  DoubleOperator     operator - (DoubleOperator lhs, DoubleOperator const& rhs);                            //!< DoubleOperator-DoubleOperator
  ///@}

  /**
   * @brief Method which prints DoubleOperator with std::cout <<.
   */
  std::ostream& operator << (std::ostream& os, DoubleOperator const& dop);
}

