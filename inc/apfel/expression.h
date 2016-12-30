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
     * @brief The default constructor, assumes that the expression does not depend on any mass index of flavour number.
     */
    Expression();

    /**
     * @brief "Massive" constructor that depends on the mass index.
     * @param MassIndex mass index.
     */
    Expression(double const& MassIndex);

    /**
     * @brief "Massless" constructor that depends on the number of active flabours.
     * @param FlavNumb number of flavours.
     */
    Expression(int const& FlavNumb);

    /**
     * @brief Virtual regular term.
     * @param x the integration variable.
     * @return the regular term at x.
     */
    virtual double Regular(double const& x) const = 0;

    /**
     * @brief Virtual singular term.
     * @param x the integration variable.
     * @return the singular term at x.
     */
    virtual double Singular(double const& x) const = 0;

    /**
     * @brief Virtual local term.
     * @param x the physical variable.
     * @return the local term at x.
     */
    virtual double Local(double const& x) const = 0;

    double MassIndex()     const { return _MassIndex; } //<! Return the mass index
    double FlavourNumber() const { return _FlavNumb; }  //<! Return the Number of flavours

  protected:
    double _MassIndex; //<! Mass index to be used in the epressions
    double _FlavNumb;  //<! Number of active flavours
  };
}
