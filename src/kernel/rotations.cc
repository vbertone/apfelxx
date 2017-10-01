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
  double PhysToQCDEv(int const& i, double const& x, double const& Q, function<double(int const&, double const&, double const&)> const& InDistFunc)
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
	+ InDistFunc(2,x,Q) + InDistFunc(-2,x,Q)
	- ( InDistFunc(1,x,Q) + InDistFunc(-1,x,Q) );
    // V3
    else if (i == 4)
      return
	+ InDistFunc(2,x,Q) - InDistFunc(-2,x,Q)
	- ( InDistFunc(1,x,Q) - InDistFunc(-1,x,Q) );
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

  //_____________________________________________________________________________
  map<int,double> PhysToQCDEv(double const& x, double const& Q, function<map<int,double>(double const&, double const&)> const& InDistFunc)
  {
    // Call function in the physical basis.
    map<int,double> PhysMap = InDistFunc(x,Q);

    // Fill in keys that don't exist.
    // Gluon (assumes that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap[0] = PhysMap[21];

    // Quarks (Fill in with zero if they don't exist).
    for (auto i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
	PhysMap[i] = 0;

    // Fill in map in the QCD evolution basis.
    map<int,double> QCDEvMap;
    QCDEvMap[0] = PhysMap.at(0);
    QCDEvMap[1] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      + PhysMap.at(5) + PhysMap.at(-5)
      + PhysMap.at(6) + PhysMap.at(-6);
    QCDEvMap[2] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      + PhysMap.at(5) - PhysMap.at(-5)
      + PhysMap.at(6) - PhysMap.at(-6);
    QCDEvMap[3] =
      + PhysMap.at(2) + PhysMap.at(-2)
      - ( PhysMap.at(1) + PhysMap.at(-1) );
    QCDEvMap[4] =
      + PhysMap.at(2) - PhysMap.at(-2)
      - ( PhysMap.at(1) - PhysMap.at(-1) );
    QCDEvMap[5] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      - 2 * ( PhysMap.at(3) + PhysMap.at(-3) );
    QCDEvMap[6] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      - 2 * ( PhysMap.at(3) - PhysMap.at(-3) );
    QCDEvMap[7] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      - 3 * ( PhysMap.at(4) + PhysMap.at(-4) );
    QCDEvMap[8] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      - 3 * ( PhysMap.at(4) - PhysMap.at(-4) );
    QCDEvMap[9] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      - 4 * ( PhysMap.at(5) + PhysMap.at(-5) );
    QCDEvMap[10] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      - 4 * ( PhysMap.at(5) - PhysMap.at(-5) );
    QCDEvMap[11] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      + PhysMap.at(5) + PhysMap.at(-5)
      - 5 * ( PhysMap.at(6) + PhysMap.at(-6) );
    QCDEvMap[12] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      + PhysMap.at(5) - PhysMap.at(-5)
      - 5 * ( PhysMap.at(6) - PhysMap.at(-6) );

    return QCDEvMap;
  }

  //_____________________________________________________________________________
  map<int,double> PhysToQCDEv(map<int,double> const& PhysMap)
  {
    // Fill in map in the QCD evolution basis. It attumes that the
    // gluon has key zero and all keys from -6 to 6 exist.
    map<int,double> QCDEvMap;
    QCDEvMap[0] = PhysMap.at(0);
    QCDEvMap[1] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      + PhysMap.at(5) + PhysMap.at(-5)
      + PhysMap.at(6) + PhysMap.at(-6);
    QCDEvMap[2] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      + PhysMap.at(5) - PhysMap.at(-5)
      + PhysMap.at(6) - PhysMap.at(-6);
    QCDEvMap[3] =
      + PhysMap.at(2) + PhysMap.at(-2)
      - ( PhysMap.at(1) + PhysMap.at(-1) );
    QCDEvMap[4] =
      + PhysMap.at(2) - PhysMap.at(-2)
      - ( PhysMap.at(1) - PhysMap.at(-1) );
    QCDEvMap[5] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      - 2 * ( PhysMap.at(3) + PhysMap.at(-3) );
    QCDEvMap[6] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      - 2 * ( PhysMap.at(3) - PhysMap.at(-3) );
    QCDEvMap[7] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      - 3 * ( PhysMap.at(4) + PhysMap.at(-4) );
    QCDEvMap[8] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      - 3 * ( PhysMap.at(4) - PhysMap.at(-4) );
    QCDEvMap[9] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      - 4 * ( PhysMap.at(5) + PhysMap.at(-5) );
    QCDEvMap[10] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      - 4 * ( PhysMap.at(5) - PhysMap.at(-5) );
    QCDEvMap[11] =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      + PhysMap.at(5) + PhysMap.at(-5)
      - 5 * ( PhysMap.at(6) + PhysMap.at(-6) );
    QCDEvMap[12] =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      + PhysMap.at(5) - PhysMap.at(-5)
      - 5 * ( PhysMap.at(6) - PhysMap.at(-6) );

    return QCDEvMap;
  }

  //_____________________________________________________________________________
  map<int,double> QCDEvToPhys(map<int,double> const& QCDEvMap)
  {
    // Rotation matrix.
    const double RotQCDEvToPhys[6][6] = {{1./12.,  1./4., 1./12., 1./24.,  1./40.,  1./60.},
					 {1./12., -1./4., 1./12., 1./24.,  1./40.,  1./60.},
					 {1./12.,     0., -1./6., 1./24.,  1./40.,  1./60.},
					 {1./12.,     0.,     0., -1./8.,  1./40.,  1./60.},
					 {1./12.,     0.,     0.,     0., -1./10.,  1./60.},
					 {1./12.,     0.,     0.,     0.,      0., -1./12.}};

    // Fill in map in the physical basis. It attumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    map<int,double> PhysMap;
    PhysMap[0] = QCDEvMap.at(0);

    // Perform the rotation.
    for (int i = 1; i <= 6; i++)
      {
	PhysMap[i]  = 0;
	PhysMap[-i] = 0;
	  for (int j = 1; j <= 6; j++)
	    {
	      PhysMap[i]  += RotQCDEvToPhys[i-1][j-1] * ( QCDEvMap.at(2*j-1) + QCDEvMap.at(2*j) );
	      PhysMap[-i] += RotQCDEvToPhys[i-1][j-1] * ( QCDEvMap.at(2*j-1) - QCDEvMap.at(2*j) );
	    }
      }

    return PhysMap;
  }
}
