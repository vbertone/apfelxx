//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/rotations.h"

using namespace std;

namespace apfel {

  //_____________________________________________________________________________
  double QCDEvToPhys(int const& i, double const& x, double const& Q, function<double(int const&, double const&, double const&)> const& InDistFunc)
  {
    // Gluon
    if      (i == 0)
      return InDistFunc(0,x,Q);
    // Singlet
    else if (i == 1)
      return
	+ InDistFunc(1,x,Q) + InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) + InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) + InDistFunc(-3,x,Q)
	+ InDistFunc(4,x,Q) + InDistFunc(-4,x,Q)
	+ InDistFunc(5,x,Q) + InDistFunc(-5,x,Q)
	+ InDistFunc(6,x,Q) + InDistFunc(-6,x,Q);
    // Valence
    else if (i == 2)
      return
	+ InDistFunc(1,x,Q) - InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) - InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) - InDistFunc(-3,x,Q)
	+ InDistFunc(4,x,Q) - InDistFunc(-4,x,Q)
	+ InDistFunc(5,x,Q) - InDistFunc(-5,x,Q)
	+ InDistFunc(6,x,Q) - InDistFunc(-6,x,Q);
    // T3
    else if (i == 3)
      return
	+ InDistFunc(1,x,Q) + InDistFunc(-1,x,Q)
	- ( InDistFunc(2,x,Q) + InDistFunc(-2,x,Q) );
    // V3
    else if (i == 4)
      return
	+ InDistFunc(1,x,Q) - InDistFunc(-1,x,Q)
	- ( InDistFunc(2,x,Q) - InDistFunc(-2,x,Q) );
    // T8
    else if (i == 5)
      return
	+ InDistFunc(1,x,Q) + InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) + InDistFunc(-2,x,Q)
	- 2 * ( InDistFunc(3,x,Q) + InDistFunc(-3,x,Q) );
    // V8
    else if (i == 6)
      return
	+ InDistFunc(1,x,Q) - InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) - InDistFunc(-2,x,Q)
	- 2 * ( InDistFunc(3,x,Q) - InDistFunc(-3,x,Q) );
    // T15
    else if (i == 7)
      return
	+ InDistFunc(1,x,Q) + InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) + InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) + InDistFunc(-3,x,Q)
	- 3 * ( InDistFunc(4,x,Q) + InDistFunc(-4,x,Q) );
    // V15
    else if (i == 8)
      return
	+ InDistFunc(1,x,Q) - InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) - InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) - InDistFunc(-3,x,Q)
	- 3 * ( InDistFunc(4,x,Q) - InDistFunc(-4,x,Q) );
    // T24
    else if (i == 9)
      return
	+ InDistFunc(1,x,Q) + InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) + InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) + InDistFunc(-3,x,Q)
	+ InDistFunc(4,x,Q) + InDistFunc(-4,x,Q)
	- 4 * ( InDistFunc(5,x,Q) + InDistFunc(-5,x,Q) );
    // V24
    else if (i == 10)
      return
	+ InDistFunc(1,x,Q) - InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) - InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) - InDistFunc(-3,x,Q)
	+ InDistFunc(4,x,Q) - InDistFunc(-4,x,Q)
	- 4 * ( InDistFunc(5,x,Q) - InDistFunc(-5,x,Q) );
    // T35
    else if (i == 11)
      return
	+ InDistFunc(1,x,Q) + InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) + InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) + InDistFunc(-3,x,Q)
	+ InDistFunc(4,x,Q) + InDistFunc(-4,x,Q)
	+ InDistFunc(5,x,Q) + InDistFunc(-5,x,Q)
	- 5 * ( InDistFunc(6,x,Q) + InDistFunc(-6,x,Q) );
    // V35
    else if (i == 12)
      return
	+ InDistFunc(1,x,Q) - InDistFunc(-1,x,Q)
	+ InDistFunc(2,x,Q) - InDistFunc(-2,x,Q)
	+ InDistFunc(3,x,Q) - InDistFunc(-3,x,Q)
	+ InDistFunc(4,x,Q) - InDistFunc(-4,x,Q)
	+ InDistFunc(5,x,Q) - InDistFunc(-5,x,Q)
	- 5 * ( InDistFunc(6,x,Q) - InDistFunc(-6,x,Q) );
    else
      return 0;
  }

}
