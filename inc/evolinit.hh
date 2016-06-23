/*
  evolinit.hh:

  Author: Valerio Bertone
*/

#pragma once

#include "evolsetup.hh"
#include "grid.hh"
#include "opset.hh"

using namespace std;

namespace apfel {

  /**
   * class evolinit:
   * class that computes the convolution of splitting functions
   * (and possibly matching conditions) and the interpolation
   * functions on the x-space interpolation grids according to
   * the settings defined in "evolsetup".
   **/
  class evolinit {
  private:
    // Attributes
    evolsetup    _Setup;
    grid         _GlobalGrid;
    opset *_sfQCD;
    //opset *_sfQED;

  public:
    // Constructors
    evolinit(evolsetup Setup_): _Setup(Setup_) { Initializer(_Setup); }
    evolinit() { Initializer(_Setup); }

    // Getters
    evolsetup Setup()      { return _Setup; }
    grid      GlobalGrid() { return _GlobalGrid; }

    // Initializer
    void Initializer(evolsetup Setup_);
  };

}
