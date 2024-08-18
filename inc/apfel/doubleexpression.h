//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @brief The DoubleExpression class encapsulates in a proper form a
   * given two-variable analytic expression in such a way that it can
   * be transformed into an double operator (only zero-mass and
   * forward expressions can be accommodated).
   */
  class DoubleExpression
  {
  public:
    virtual ~DoubleExpression() = default;

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The "DoubleExpression" constructor
     * @param expr1: first expression
     * @param expr2: second expression
     */
    DoubleExpression(Expression const& expr1 = Expression{}, Expression const& expr2 = Expression{});
    ///@}

    /**
     * @name Expression components The different possible components
     * of an expression: local-local, local-singular, local-regular,
     * singular-local, singular-singular, singular-regular,
     * regular-local, regular-singular, regular-regular.
     */
    ///@{
    /**
     * @brief Virtual local-local term.
     * @return The local-local term at x
     */
    virtual double LocalLocal(double const& x1, double const& x2) const { return _expr1.Local(x1) * _expr2.Local(x2); }

    /**
     * @brief Virtual local-singular term.
     * @return The local-singular term at x
     */
    virtual double LocalSingular(double const& x1, double const& x2) const { return _expr1.Local(x1) * _expr2.Singular(x2); }

    /**
     * @brief Virtual local-regular term.
     * @return The local-regular term at x
     */
    virtual double LocalRegular(double const& x1, double const& x2) const { return _expr1.Local(x1) * _expr2.Regular(x2); }

    /**
     * @brief Virtual singular-local term.
     * @return The singular-local term at x
     */
    virtual double SingularLocal(double const& x1, double const& x2) const { return _expr1.Singular(x1) * _expr2.Local(x2); }

    /**
     * @brief Virtual singular-singular term.
     * @return The singular-singular term at x
     */
    virtual double SingularSingular(double const& x1, double const& x2) const { return _expr1.Singular(x1) * _expr2.Singular(x2); }

    /**
     * @brief Virtual singular-regular term.
     * @return The singular-regular term at x
     */
    virtual double SingularRegular(double const& x1, double const& x2) const { return _expr1.Singular(x1) * _expr2.Regular(x2); }

    /**
     * @brief Virtual regular-local term.
     * @return The regular-local term at x
     */
    virtual double RegularLocal(double const& x1, double const& x2) const { return _expr1.Regular(x1) * _expr2.Local(x2); }

    /**
     * @brief Virtual regular-singular term.
     * @return The regular-singular term at x
     */
    virtual double RegularSingular(double const& x1, double const& x2) const { return _expr1.Regular(x1) * _expr2.Singular(x2); }

    /**
     * @brief Virtual regular-regular term.
     * @return The regular-regular term at x
     */
    virtual double RegularRegular(double const& x1, double const& x2) const { return _expr1.Regular(x1) * _expr2.Regular(x2); }
    ///@}

  private:
    Expression const& _expr1; //!< First expression
    Expression const& _expr2; //!< Second expression
  };
  ///@}
}
