//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @brief The DoubleExpression class encapsulates in a proper form a
   *  two-variable analytic expression in such a way that it can be
   *  transformed into an DoubleOperator object. This is meant to be
   *  used with the partonic cross sections of processes such as SIDIS
   *  and Drell-Yan that require double Mellin convolutions.
   * @note For now only zero-mass and forward expressions can be
   *  accommodated.
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
     */
    DoubleExpression();
    ///@}

    /**
     * @name DoubleExpression components
     * The different possible components of a double expression.
     */
    ///@{
    /**
     * @brief Virtual function for the local-local term.
     * @return The local-local term at x1 and x2
     */
    virtual double LocalLocal(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the local-singular term.
     * @return The local-singular term at x1 and x2
     */
    virtual double LocalSingular(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the local-regular term.
     * @return The local-regular term at x1 and x2
     */
    virtual double LocalRegular(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the singular-local term.
     * @return The singular-local term at x1 and x2
     */
    virtual double SingularLocal(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the singular-singular term.
     * @return The singular-singular term at x1 and x2
     */
    virtual double SingularSingular(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the singular-regular term.
     * @return The singular-regular term at x1 and x2
     */
    virtual double SingularRegular(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the regular-local term.
     * @return The regular-local term at x1 and x2
     */
    virtual double RegularLocal(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the regular-singular term.
     * @return The regular-singular term at x1 and x2
     */
    virtual double RegularSingular(double const&, double const&) const { return 0; }

    /**
     * @brief Virtual function for the regular-regular term.
     * @return The regular-regular term at x1 and x2
     */
    virtual double RegularRegular(double const&, double const&) const { return 0; }
    ///@}
  };
}
