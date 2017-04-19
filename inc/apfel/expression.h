//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

namespace apfel
{
  /**
   * @brief The Expression class for the manipulation of the splitting and coeffient functions.
   *
   * This class encapsulate in a proper form an analystic expressions in such a way that
   * it can be transformed into an operator.
   */
  class Expression
  {
  public:
    /**
     * @brief The default constructor
     */
    Expression();

    /**
     * @brief Virtual regular term.
     * @param x the integration variable.
     * @return the regular term at x.
     */
    virtual double Regular(double const&) const { return 0; }

    /**
     * @brief Virtual singular term.
     * @param x the integration variable.
     * @return the singular term at x.
     */
    virtual double Singular(double const&) const { return 0; }

    /**
     * @brief Virtual local term.
     * @param x the physical variable.
     * @return the local term at x.
     */
    virtual double Local(double const&) const { return 0; }
  };

  /**
   * @brief Identity expression (delta function)
   */
  class Identity: public Expression
  {
  public:
  Identity(): Expression() { }
    double Local(double const&) const { return 1; }
  };

  /**
   * @brief Zero expression
   */
  class Null: public Expression
  {
  public:
  Null(): Expression() { }
  };

}
