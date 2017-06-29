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
  unordered_map<int,double> PhysToQCDEv(double const& x, double const& Q, function<map<int,double>(double const&, double const&)> const& InDistFunc)
  {
    // Call function in the physical basis.
    map<int,double> PhysMap = InDistFunc(x,Q);

    // Fill in keys that don't exist.
    // Gluon (assumes that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap.insert({0,PhysMap.at(21)});

    // Quarks (Fill in with zero if they don't exist).
    for (auto i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
	PhysMap.insert({i,0});

    // Define distributions in the QCD evolution basis.
    const double Gluon = PhysMap.at(0);
    const double Singlet =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      + PhysMap.at(5) + PhysMap.at(-5)
      + PhysMap.at(6) + PhysMap.at(-6);
    const double Valence =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      + PhysMap.at(5) - PhysMap.at(-5)
      + PhysMap.at(6) - PhysMap.at(-6);
    const double T3 =
      + PhysMap.at(2) + PhysMap.at(-2)
      - ( PhysMap.at(1) + PhysMap.at(-1) );
    const double V3 =
      + PhysMap.at(2) - PhysMap.at(-2)
      - ( PhysMap.at(1) - PhysMap.at(-1) );
    const double T8 =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      - 2 * ( PhysMap.at(3) + PhysMap.at(-3) );
    const double V8 =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      - 2 * ( PhysMap.at(3) - PhysMap.at(-3) );
    const double T15 =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      - 3 * ( PhysMap.at(4) + PhysMap.at(-4) );
    const double V15 =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      - 3 * ( PhysMap.at(4) - PhysMap.at(-4) );
    const double T24 =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      - 4 * ( PhysMap.at(5) + PhysMap.at(-5) );
    const double V24 =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      - 4 * ( PhysMap.at(5) - PhysMap.at(-5) );
    const double T35 =
      + PhysMap.at(1) + PhysMap.at(-1)
      + PhysMap.at(2) + PhysMap.at(-2)
      + PhysMap.at(3) + PhysMap.at(-3)
      + PhysMap.at(4) + PhysMap.at(-4)
      + PhysMap.at(5) + PhysMap.at(-5)
      - 5 * ( PhysMap.at(6) + PhysMap.at(-6) );
    const double V35 =
      + PhysMap.at(1) - PhysMap.at(-1)
      + PhysMap.at(2) - PhysMap.at(-2)
      + PhysMap.at(3) - PhysMap.at(-3)
      + PhysMap.at(4) - PhysMap.at(-4)
      + PhysMap.at(5) - PhysMap.at(-5)
      - 5 * ( PhysMap.at(6) - PhysMap.at(-6) );

    // Fill in map in the QCD evolution basis.
    unordered_map<int,double> QCDEvMap;
    QCDEvMap.insert({0 , Gluon});
    QCDEvMap.insert({1 , Singlet});
    QCDEvMap.insert({2 , Valence});
    QCDEvMap.insert({3 , T3});
    QCDEvMap.insert({4 , V3});
    QCDEvMap.insert({5 , T8});
    QCDEvMap.insert({6 , V8});
    QCDEvMap.insert({7 , T15});
    QCDEvMap.insert({8 , V15});
    QCDEvMap.insert({9 , T24});
    QCDEvMap.insert({10, V24});
    QCDEvMap.insert({11, T35});
    QCDEvMap.insert({12, V35});

    return QCDEvMap;
  }

}
