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
   * analytic expression in such a way that it can be transformed into
   * an operator.
   */
  class Expression
  {
  public:
    virtual ~Expression() = default;

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The "Expression" constructor
     * @param eta: upper limit of the convolution integral (default: 1)
     */
    Expression(double const& eta = 1);
    ///@}

    /**
     * @name Expression components
     * The different possible components of an expression: regular,
     * singular, and local.
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

    /**
     * @brief Virtual local term for principal-valued integrals a la
     * ERBL with singularity at x = 1, i.e. corresponding to the
     * ++-prescription.
     * @return The local term for ++-prespribed distributions at x
     */
    virtual double LocalPP(double const&) const { return 0; }

    /**
     * @brief Virtual singular term for principal-valued integrals in
     * the DGLAP region (i.e. with pole in x in the interval (0,1)).
     * @return The singular term for principal-valued distributions at x
     */
    virtual double SingularPV(double const&) const { return 0; }

    /**
     * @brief Virtual local term for principal-valued integrals a la
     * DGLAP with singularity in the interval (0,1).
     * @return The local term for principal-valued distributions at x
     */
    virtual double LocalPV(double const&) const { return 0; }
    ///@}

    /**
     * @brief Function that sets the value of a possible external
     * variable.
     */
    void SetExternalVariable(double const& extvar) const { _extvar = extvar; }

    /**
     * @brief Function that returns the value of the scaling parameter
     * eta.
     */
    double eta() const { return _eta; }

  protected:
    double mutable _extvar;  //!< External kinematic variable
    double const   _eta;     //!< Scaling parameter
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
