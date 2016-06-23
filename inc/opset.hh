/*
  opset.hh:

  Author: Valerio Bertone
*/

#pragma once

#include <string>

#include "grid.hh"

using namespace std;

namespace apfel {

  // Kinds of contribution to the functions
  enum behaviour {
    REGULAR,
    SINGULAR,
    LOCAL
  };

  /**
   * class opset:
   * Parent class for the management of the convolutions.
   **/
  class opset {
  private:
    // Attributes
    int               _GaussPoints;
    grid              _GlobalGrid;
    vector<double***> _Operators;

  public:
    // Standard constructor
    opset(int  const& GaussPoints_,
	  grid const& GlobalGrid_):
      _GaussPoints(GaussPoints_),
      _GlobalGrid(GlobalGrid_)
    { }

    // Multiplication constructor
    opset(int const& im1, opset const* op1, int const& im2, opset const* op2);

    // Constructor to join two operators
    opset(opset const* op1, opset const* op2);

    // Destructor
    ~opset();

    // Getters
    int  GaussPoints()            const { return _GaussPoints; }
    grid GlobalGrid()             const { return _GlobalGrid; }
    vector<double***> Operators() const { return _Operators;}

    // Virtual getters
    virtual int NumberOfMembers() const = 0;

    // Virtual function to be overriden by the daughter class
    virtual double OperatorFunction(int const& member_, behaviour const& behaviour_, int const& nf_, double const& x_) const = 0;

    // Virtual map between operators and distributions
    virtual vector<int> Map(int const& id) const { return {}; }

    // Function that translates the operators on the gird
    void CreateOperators(int const& imass);

  };

}
