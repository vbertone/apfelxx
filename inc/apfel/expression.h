//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @brief The Expression class encapsulates in a proper form a given
   * analystic expression in such a way that it can be transformed
   * into an operator.
   */
  class Expression
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The "Expression" constructor
     * @param ext: switch to know if depends explicitly on the external variable (default: false)
     * @param eta: upper limit of the convolution integral (default: 1)
     */
    Expression(bool const& ext = false, double const& eta = 1);

    /**
     * @brief The "Expression" constructor
     * @param eta: upper limit of the convolution integral
     */
    Expression(double const& eta);
    ///@}

    /**
     * @name Expression components
     * The three different possible components of an expression:
     * regular, singular, and local.
     */
    ///@{
    /**
     * @brief Virtual regular term.
     * @return The regular term at x
     */
    virtual double Regular(double const&) const { return 0; }

    /**
     * @brief Virtual singular term.
     * @return The singular term at x
     */
    virtual double Singular(double const&) const { return 0; }

    /**
     * @brief Virtual local term.
     * @return The local term at x
     */
    virtual double Local(double const&) const { return 0; }
    ///@}

    /**
     * @brief Function that returns the value of the scaling parameter
     * eta.
     */
    double eta() const { return _eta; }

    /**
     * @brief Function that returns the switch that encodes if the
     * expression depends explicitly on the external variable.
     */
    bool DependsOnExternalVariable() const { return _ext; }

    /**
     * @brief Function that sets the value of a possible external
     * variable.
     */
    void SetExternalVariable(double const& extvar) const { _extvar = extvar; }

  protected:
    bool   const   _ext;     //!< If the expression depends explicitly on the external variable
    double mutable _extvar;  //!< External kinematic variable
    double const   _eta;     //!< Mass parameter
  };

  /**
   * @defgroup RecExprs Recurrent expressions
   * Collection of recurrent expressions. This includes the identity
   * and the null expressions.
   */
  ///@{
  /**
   * @brief Derived class from Expression to implement the Identity
   * operator (delta function).
   */
  class Identity: public Expression
  {
  public:
    Identity(): Expression() { }
    double Local(double const&) const { return 1; }
  };

  /**
   * @brief Derived class from Expression to implement the Null
   * operator (zero).
   */
  class Null: public Expression
  {
  public:
    Null(): Expression() { }
  };
  ///@}
}
