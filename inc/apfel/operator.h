//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/grid.h>
#include <apfel/expression.h>
#include <apfel/distribution.h>
#include <apfel/integrator.h>
#include <apfel/lagrangeinterpolator.h>
#include <apfel/matrix.h>

namespace apfel
{
  /**
   * @brief The Operator class.
   *
   * This class defines the basic object "Operator" which is
   * essentially the convolution on the grid bewteen a function
   * (e.g. a splitting function) and the inetrpolant functions.
   */
  class Operator: protected Integrator, protected LagrangeInterpolator
  {
  public:

    Operator() = delete;

    /**
     * @brief The default constructor, takes a Grid and a vector of pointers to functions.
     * @param gr "Grid" object
     * @param func vector of pointers to one-dimentional functions
     */
    Operator(Grid const& gr, Expression const& expr, double const& eps = 1e-5);

    // operators
    Distribution operator *= (Distribution const& d) const; //!< this *= Distribution
    Operator&    operator  = (Operator const& o);           //!< this  = Operator
    Operator&    operator *= (Operator const& o);           //!< this *= Operator
    Operator&    operator *= (double const& s);             //!< this *= Scalar
    Operator&    operator /= (double const& s);             //!< this /= Scalar
    Operator&    operator += (Operator const& o);           //!< this += Operator
    Operator&    operator -= (Operator const& o);           //!< this -= Operator

  protected:
    /**
     * @see Integrator::integrand
     */
    double integrand(double const& x) const;

  private:
    Grid                   const& _grid;         //!< Grid on which to compute the operator
    Expression const*      const  _expr;         //!< Expression to be commuted into an operator
    double                        _eps;          //!< Precision of the dgauss integration
    vector<matrix<double>>        _Operator;     //!< Operator values.

    // Global variables
    int    _alpha, _beta;
    int    _ig;
    double _ws;
  };

  // Extra operation definitions where Operator is at the left hand
  // side (lhs).
  Distribution operator * (Operator lhs, Distribution const& rhs); //!< Operator*Distribution
  Operator     operator * (Operator lhs, Operator const& rhs);     //!< Operator*Operator
  Operator     operator * (double const& s, Operator rhs);         //!< Scalar*Operator
  Operator     operator * (Operator lhs, double const& s);         //!< Operator*Scalar
  Operator     operator / (Operator lhs, double const& s);         //!< Operator/Scalar
  Operator     operator + (Operator lhs, Operator const& rhs);     //!< Operator+Operator
  Operator     operator - (Operator lhs, Operator const& rhs);     //!< Operator-Operator
}
