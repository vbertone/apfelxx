/*
  splitting.hh:

  Author: Valerio Bertone
*/

#pragma once

#include "opset.hh"
#include "grid.hh"

using namespace std;

namespace apfel {

  // ======================================================================

  /**
   * class QCD_space_unpol:
   * Daugther class for the QCD unpolarized space-like splitting functions.
   **/
  class QCD_space_unpol: public opset {
  private:
    // Attributes
    string _Name;
    int    _NumberOfMembers;
    int    _GaussPoints;
    grid   _GlobalGrid;

    // Internal member functions
    // LO splitting functions
    double P0ns(behaviour const& behaviour_, int const& nf_, double const& x_) const;
    double P0qg(behaviour const& behaviour_, int const& nf_, double const& x_) const;
    double P0gq(behaviour const& behaviour_, int const& nf_, double const& x_) const;
    double P0gg(behaviour const& behaviour_, int const& nf_, double const& x_) const;

  public:
    // Constructor
    QCD_space_unpol(int const& GaussPoints_, grid const& GlobalGrid_):
      _Name("QCD space-like unpolarized splitting functions"),
      _NumberOfMembers(4),
      _GaussPoints(GaussPoints_),
      _GlobalGrid(GlobalGrid_),
      opset(GaussPoints_, GlobalGrid_)
    { }

    // Getters
    int    NumberOfMembers() const { return _NumberOfMembers; };

    // Splitting function
    double OperatorFunction(int const& member_, behaviour const& behaviour_, int const& nf_, double const& x) const;

    // Map to define how the operator members have to be combined with a set
    // of distibutions (in the evolution basis).
    vector<int> Map(int const& id) const;

  };

  // ======================================================================

  /**
   * class QCD_space_pol:
   * Daugther class for the QCD polarized space-like splitting functions.
   **/

  // ======================================================================

  /**
   * class QCD_time_unpol:
   * Derived class for the QCD unpolarized time-like splitting functions.
   **/

}
