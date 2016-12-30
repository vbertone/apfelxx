//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <apfel/integrator.h>
#include <apfel/lagrangeinterpolator.h>

namespace apfel
{
  class Grid;
  class Expression;
  template<class T> using vector3d = vector<vector<vector<T>>>;

  /**
   * @brief The Operator class.
   *
   * This class defines the basic object "Operator" which is essentially
   * the convolution on the grid bewteen a function (i.e. a splitting function)
   * and the inetrpolant functions.
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

  protected:
    /**
     * @see Integrator::integrand
     */
    double integrand(double const& x) const;

  private:    
    Grid         const& _grid;         //!< Grid on which to compute the operator
    Expression   const& _expr;         //!< Expression to be commuted into an operator
    double       const& _eps;          //!< precision of the dgauss integration
    vector3d<double> _Operator;        //!< the operator values.

    // Global variables
    int     _alpha, _beta;
    int     _ig;
    double  _ws;

  };
}
