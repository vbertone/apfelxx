//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"
#include "apfel/distribution.h"
#include "apfel/matrix.h"

namespace apfel
{
  /**
   * @brief This class defines the basic object Operator which
   * essentially contains the convolution bewteen an "Expression"
   * object (e.g. a splitting function) and a set of interpolanting
   * functions defined on a grid.
   * @note Both the forward and the non-forward cases can be handled
   * by this class. The correspond respectively to inclusive
   * (DGLAP-like) and exclusive (GPD-like) processes.
   */
  class Operator
  {
  public:
    Operator() = delete;
    Operator(Operator const&) = default;

    /**
     * @brief The Operator constructor.
     * @param gr: the Grid object
     * @param expr: the expression to be convoluted
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     * @param gpd: whether the operator had to computed for a GPD-like expression (default: false)
     */
    Operator(Grid const& gr, Expression const& expr, double const& eps = 1e-5, bool const& gpd = false);

    /**
     * @brief The Operator virtual destructor.
     */
    virtual ~Operator() {}

    /**
     * @name Binary operators
     */
    ///@{
    Distribution operator *= (Distribution const& d) const;         //!< this *= Distribution
    Operator& operator *= (Operator const& o);                      //!< this *= Operator
    Operator& operator  = (Operator const& o);                      //!< this  = Operator
    Operator& operator *= (double const& s);                        //!< this *= Scalar
    Operator& operator *= (std::function<double(double const&)> f); //!< This *= Function
    Operator& operator /= (double const& s);                        //!< this /= Scalar
    Operator& operator += (Operator const& o);                      //!< this += Operator
    Operator& operator -= (Operator const& o);                      //!< this -= Operator
    ///@}

    /**
     * @brief Function that builds a DGLAP-like operator.
     */
    void BuildOperatorDGLAP();

    /**
     * @brief Function that builds a GPD-like operator.
     */
    void BuildOperatorGPD();

    /**
     * @brief Function that builds a DGLAP-like operator, that is xdependent.
     */
    void BuildOperatorDGLAPxdep();

    /**
     * @brief Function that interpolates the operator over the first
     * index return a Distribution object.
     * @param x: the value in x to be interpolated
     */
    Distribution Evaluate(double const& x) const;

    /**
     * @brief Function that returns the Grid object associated with
     * the operator.
     */
    Grid const& GetGrid() const { return _grid; }

    /**
     * @brief Function that returns the Expression object associated
     * with the operator.
     */
    Expression const& GetExpression() const { return _expr; }

    /**
     * @brief Function that returns the integration accuracy,
     */
    double const& GetIntegrationAccuracy() const { return _eps; }

    /**
     * @brief Function that returns the GPD switch,
     */
    bool const& IsGPD() const { return _gpd; }

    /**
     * @brief Function that returns the Operator container.
     */
    std::vector<matrix<double>> GetOperator() const { return _Operator; }

    /**
     * @brief Function that prints the Operator object.
     */
    void Print() const { std::cout << *this << std::endl; }

  protected:
    Grid                 const& _grid;      //!< Grid on which to compute the operator
    Expression           const& _expr;      //!< Expression to be used
    double               const  _eps;       //!< Integration accuracy
    bool                 const  _gpd;       //!< GPD switch
    std::vector<matrix<double>> _Operator;  //!< Operator values

    friend std::ostream& operator << (std::ostream& os, Operator const& op);
  };

  /**
   * @name Ternary operators
   */
  ///@{
  Distribution operator * (Operator lhs, Distribution const& rhs);                //!< Operator*Distribution
  Operator     operator * (Operator lhs, Operator const& rhs);                    //!< Operator*Operator
  Operator     operator * (double const& s, Operator rhs);                        //!< Scalar*Operator
  Operator     operator * (Operator lhs, double const& s);                        //!< Operator*Scalar
  Operator     operator * (std::function<double(double const&)> f, Operator rhs); //!< function*Operator
  Operator     operator * (Operator lhs, std::function<double(double const&)> f); //!< Operator*function
  Operator     operator / (Operator lhs, double const& s);                        //!< Operator/Scalar
  Operator     operator + (Operator lhs, Operator const& rhs);                    //!< Operator+Operator
  Operator     operator - (Operator lhs, Operator const& rhs);                    //!< Operator-Operator
  ///@}

  /**
   * @brief Method which prints Operator with std::cout <<.
   */
  std::ostream& operator << (std::ostream& os, Operator const& op);
}
