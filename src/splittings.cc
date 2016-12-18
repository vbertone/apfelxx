/*
  splittings.cc:

  Author: Valerio Bertone
*/

#include <iostream>
#include <cmath>

#include "apfel/colorfactors.h"
#include "apfel/splittings.h"

using namespace std;

namespace apfel {

  /**
   * Unpolarized space-like splitting functions.
   **/
  // ================================================================================
  // Splitting function
  double QCD_space_unpol::OperatorFunction(int const& member, behaviour const& behaviour, int const& nf, double const& x) const
  {
    double sf = 0;

    // Perturbative order
    switch(member)
      {
	// LO
      case(0):
	sf = P0ns(behaviour,nf,x);
	break;
      case(1):
	sf = P0qg(behaviour,nf,x);
	break;
      case(2):
	sf = P0qg(behaviour,nf,x);
	break;
      case(3):
	sf = P0gg(behaviour,nf,x);
	break;
      default:
	sf = 0;
	break;
      }
    return sf;
  };

  // ================================================================================
  // LO splitting functions
  double QCD_space_unpol::P0ns(behaviour const& behaviour, int const& nf, double const& x) const
  {
    switch(behaviour)
      {
      case(REGULAR):
	return - 2 * CF * ( 1 + x );
	break;
      case(SINGULAR):
	return 4 * CF / ( 1 - x );
	break;
      case(LOCAL):
	return 4 * CF * log( 1 - x ) + 3 * CF;
	break;
      default:
	return 0;
	break;
      }
  };

  double QCD_space_unpol::P0qg(behaviour const& behaviour, int const& nf, double const& x) const
  {
    switch(behaviour)
      {
      case(REGULAR):
	return 2 * nf * ( 1 - 2 * x + 2 * x * x );
	break;
      default:
	return 0;
	break;
      }
  };

  double QCD_space_unpol::P0gq(behaviour const& behaviour, int const& nf, double const& x) const
  {
    switch(behaviour)
      {
      case(REGULAR):
	return 4 * CF * ( - 1 + 0.5 * x + 1 / x );
	break;
      default:
	return 0;
	break;
      }
  };

  double QCD_space_unpol::P0gg(behaviour const& behaviour, int const& nf, double const& x) const
  {
    switch(behaviour)
      {
      case(REGULAR):
	return 4 * CA * ( - 2 + x - x * x + 1 / x );
	break;
      case(SINGULAR):
	return 4 * CA / ( 1 - x );
	break;
      case(LOCAL):
	return 4 * CA * log( 1 - x ) - 2 / 3 * nf + 11 / 3 * CA;
	break;
      default:
	return 0;
	break;
      }
  };

  // Map to define how the operator members have to be combined with a set
  // of distibutions (in the evolution basis).
  vector<int> QCD_space_unpol::Map(int const& id) const
  {
    if     (id == 0)  return {};     // Photon
    else if(id == 1)  return {0, 1}; // Singlet
    else if(id == 2)  return {2, 3}; // Gluon
    else if(id == 3)  return {0};    // T3
    else if(id == 4)  return {0};    // T8
    else if(id == 5)  return {0};    // T15
    else if(id == 6)  return {0};    // T24
    else if(id == 7)  return {0};    // T35
    else if(id == 8)  return {0};    // Total valence
    else if(id == 9)  return {0};    // V3
    else if(id == 10) return {0};    // V8
    else if(id == 11) return {0};    // V15
    else if(id == 12) return {0};    // V24
    else if(id == 13) return {0};    // V35
    else              return {};
  }

}
