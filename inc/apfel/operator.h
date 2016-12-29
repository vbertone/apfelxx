//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include <vector>
#include <memory>

#include "apfel/integrator.h"

using std::vector;
using std::unique_ptr;

namespace apfel
{

  class SubGrid;
  class Grid;
  class Expression;
  class Interpolator;

  /**
   * @brief The Operator class.
   *
   * This class defines the basic object "Operator" which is essentially
   * the convolution on the grid bewteen a function (i.e. a splitting function)
   * and the inetrpolant functions.
   */
  class Operator: public Integrator
  {
  public:

    Operator() = delete;

    /**
     * @brief The default constructor, takes a Grid and a vector of pointers to functions.
     * @param gr "Grid" object
     * @param func vector of pointers to one-dimentional functions
     */
    Operator(Grid const& gr, Expression const& expr, Interpolator const& inter, double const& eps = 1e-5);  //!< The class constructor

    /**
     * @see Integrator::integrand
     */
    double integrand(double const& x) const;

  private:
    Grid         const& _grid;         //!< Grid on which to compute the operator
    Expression   const& _expr;         //!< Expression to be commuted into an operator
    Interpolator const& _inter;        //!< Object cobstaining the interpolants
    double       const& _eps;          //!< precision of the dgauss integration
    double       ***_Operator;

    // Global variables
    int     _alpha, _beta;
    int     _ig;
    double  _ws;

  };
}
