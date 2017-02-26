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
    virtual double Regular(double const& x) const { return 0 * x; };

    /**
     * @brief Virtual singular term.
     * @param x the integration variable.
     * @return the singular term at x.
     */
    virtual double Singular(double const& x) const { return 0 * x; };

    /**
     * @brief Virtual local term.
     * @param x the physical variable.
     * @return the local term at x.
     */
    virtual double Local(double const& x) const { return 0 * x; };
  };
}
