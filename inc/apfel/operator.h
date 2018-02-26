//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/grid.h"
#include "apfel/expression.h"
#include "apfel/distribution.h"
#include "apfel/integrator.h"
#include "apfel/lagrangeinterpolator.h"
#include "apfel/matrix.h"

typedef double accuracy;

namespace apfel
{
  /**
   * @brief The Operator class defines the basic object "Operator"
   * which is essentially the convolution on the grid bewteen an
   * Expression object (e.g. a splitting function) and the inetrpolant
   * functions.
   */
  class Operator: protected Integrator, protected LagrangeInterpolator
  {
  public:

    Operator() = delete;

    /**
     * @brief The Operator default constructor.
     * @param gr: the Grid object
     * @param expr: the expression to be transformed
     * @param eps: relative accuracy of the numerical integrations (default: 10<SUP>-5</SUP>)
     */
    Operator(Grid const& gr, Expression const& expr, accuracy const& eps = 1e-5);

    /**
     * @name Operator binary operators
     * Binary operators involving an Operator.
     */
    ///@{
    Distribution operator *= (Distribution const& d) const;       //!< this *= Distribution
    Operator&    operator  = (Operator const& o);                 //!< this  = Operator
    Operator&    operator *= (Operator const& o);                 //!< this *= Operator
    Operator&    operator *= (double const& s);                   //!< this *= Scalar
    Operator&    operator *= (function<double(double const&)> f); //!< This *= function
    Operator&    operator /= (double const& s);                   //!< this /= Scalar
    Operator&    operator += (Operator const& o);                 //!< this += Operator
    Operator&    operator -= (Operator const& o);                 //!< this -= Operator
    ///@}

    /**
     * @brief Function that returns the Grid object of the operator.
     */
    Grid const& GetGrid() const { return _grid; }

  protected:
    /**
     * @see Integrator::integrand
     */
    double integrand(double const& x) const;

  private:
    Grid                   const& _grid;         //!< Grid on which to compute the operator
    Expression const*      const  _expr;         //!< Expression to be commuted into an operator
    vector<matrix<double>>        _Operator;     //!< Operator values.

    // Global variables
    int    _alpha;
    int    _ig;
    double _ws;
    double _xbeta;
    double _eta;
  };

  /**
   * @name Operator ternary operators
   * Ternary operators involving Operators.
   */
  ///@{
  Distribution operator * (Operator lhs, Distribution const& rhs);           //!< Operator*Distribution
  Operator     operator * (Operator lhs, Operator const& rhs);               //!< Operator*Operator
  Operator     operator * (double const& s, Operator rhs);                   //!< Scalar*Operator
  Operator     operator * (Operator lhs, double const& s);                   //!< Operator*Scalar
  Operator     operator * (function<double(double const&)> f, Operator rhs); //!< function*Operator
  Operator     operator * (Operator lhs, function<double(double const&)> f); //!< Operator*function
  Operator     operator / (Operator lhs, double const& s);                   //!< Operator/Scalar
  Operator     operator + (Operator lhs, Operator const& rhs);               //!< Operator+Operator
  Operator     operator - (Operator lhs, Operator const& rhs);               //!< Operator-Operator
  ///@}
}
