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
  class Distribution;
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

    /**
     * @brief Copy constructor with precomputed operator
     * @param obj the input copy object
     * @param op the precomputed operator
     */
    Operator(Operator const& obj, vector3d<double> const& op);

    // operators
    Distribution operator*(Distribution const& d) const; //!< Operator * Distribution
    Operator     operator*(Operator const& o)     const; //!< Operator * Operator

  protected:
    /**
     * @see Integrator::integrand
     */
    double integrand(double const& x) const;

  private:    
    Grid             const& _grid;         //!< Grid on which to compute the operator
    Expression       const& _expr;         //!< Expression to be commuted into an operator
    double           const& _eps;          //!< Precision of the dgauss integration
    vector3d<double>        _Operator;     //!< Operator values.

    // Global variables
    int    _alpha, _beta;
    int    _ig;
    double _ws;
  };
}
