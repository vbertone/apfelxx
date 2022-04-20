//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/rotations.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, double> PhysToQCDEv(std::map<int, double> const& InPhysMap)
  {
    // Call function in the physical basis.
    std::map<int, double> PhysMap = InPhysMap;

    // Fill in keys that don't exist.
    // Gluon (assumes that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap[0] = PhysMap[21];

    // Quarks (Fill in with zero if they don't exist).
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap[i] = 0;

    // Fill in map in the QCD evolution basis. It attumes that the
    // gluon has key zero and all keys from -6 to 6 exist.
    std::map<int,double> QCDEvMap;
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
  std::map<int, double> QCDEvToPhys(std::map<int, double> const& QCDEvMap)
  {
    // Fill in map in the physical basis. It attumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, double> PhysMap;
    PhysMap[0]  = QCDEvMap.at(0);
    PhysMap[21] = QCDEvMap.at(0);

    // Perform the rotation.
    for (int i = 1; i <= 6; i++)
      {
        PhysMap[i]  = 0;
        PhysMap[-i] = 0;
        for (int j = 1; j <= 6; j++)
          {
            PhysMap[i]  += ( RotQCDEvToPhys[i-1][j-1] / 2 ) * ( QCDEvMap.at(2*j-1) + QCDEvMap.at(2*j) );
            PhysMap[-i] += ( RotQCDEvToPhys[i-1][j-1] / 2 ) * ( QCDEvMap.at(2*j-1) - QCDEvMap.at(2*j) );
          }
      }
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, Distribution> QCDEvToPhys(std::map<int, Distribution> const& QCDEvMap)
  {
    // Fill in map in the physical basis. It attumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, Distribution> PhysMap;
    PhysMap.insert({0, QCDEvMap.at(0)});

    // Perform the rotation.
    for (int i = 1; i <= 6; i++)
      {
        Distribution Td = ( RotQCDEvToPhys[i-1][0] / 2 ) * ( QCDEvMap.at(1) + QCDEvMap.at(2) );
        Distribution Vd = ( RotQCDEvToPhys[i-1][0] / 2 ) * ( QCDEvMap.at(1) - QCDEvMap.at(2) );
        for (int j = 2; j <= 6; j++)
          {
            Td += ( RotQCDEvToPhys[i-1][j-1] / 2 ) * ( QCDEvMap.at(2*j-1) + QCDEvMap.at(2*j) );
            Vd += ( RotQCDEvToPhys[i-1][j-1] / 2 ) * ( QCDEvMap.at(2*j-1) - QCDEvMap.at(2*j) );
          }
        PhysMap.insert({i, Td});
        PhysMap.insert({-i, Vd});
      }
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, Operator> QCDEvToPhys(std::map<int, Operator> const& QCDEvMap)
  {
    // Fill in map in the physical basis. It attumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, Operator> PhysMap;
    PhysMap.insert({0, QCDEvMap.at(0)});

    // Perform the rotation.
    for (int i = 1; i <= 6; i++)
      {
        Operator Td = ( RotQCDEvToPhys[i-1][0] / 2 ) * ( QCDEvMap.at(1) + QCDEvMap.at(2) );
        Operator Vd = ( RotQCDEvToPhys[i-1][0] / 2 ) * ( QCDEvMap.at(1) - QCDEvMap.at(2) );
        for (int j = 2; j <= 6; j++)
          {
            Td += ( RotQCDEvToPhys[i-1][j-1] / 2 ) * ( QCDEvMap.at(2*j-1) + QCDEvMap.at(2*j) );
            Vd += ( RotQCDEvToPhys[i-1][j-1] / 2 ) * ( QCDEvMap.at(2*j-1) - QCDEvMap.at(2*j) );
          }
        PhysMap.insert({i, Td});
        PhysMap.insert({-i, Vd});
      }
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, double> PhysToPlusMinus(std::map<int, double> const& InPhysMap)
  {
    // Call function in the physical basis.
    std::map<int, double> PhysMap = InPhysMap;

    // Fill in keys that don't exist.
    // Gluon (assumes that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap[0] = PhysMap[21];

    // Quarks (Fill in with zero if they don't exist).
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap[i] = 0;

    // Fill in map in the PlusMinus basis. It attumes that the gluon
    // has key zero and all keys from -6 to 6 exist.
    std::map<int,double> QCDEvMap;
    QCDEvMap[0]  = PhysMap.at(0);
    for (int i = 1; i <= 6; i++)
      {
        QCDEvMap[2*i-1] = PhysMap.at(i) + PhysMap.at(-i);
        QCDEvMap[2*i]   = PhysMap.at(i) - PhysMap.at(-i);
      }

    return QCDEvMap;
  }

  //_____________________________________________________________________________
  std::map<int, double> PlusMinusToPhys(std::map<int, double> const& PlusMinusMap)
  {
    // Fill in map in the physical basis. It attumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, double> PhysMap;
    PhysMap[0]  = PlusMinusMap.at(0);
    PhysMap[21] = PlusMinusMap.at(0);

    // Fill in map in the physical basis.
    for (int i = 1; i <= 6; i++)
      {
        PhysMap[i]  = ( PlusMinusMap.at(2*i-1) + PlusMinusMap.at(2*i) ) / 2;
        PhysMap[-i] = ( PlusMinusMap.at(2*i-1) - PlusMinusMap.at(2*i) ) / 2;
      }

    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, Distribution> PlusMinusToPhys(std::map<int, Distribution> const& PlusMinusMap)
  {
    // Fill in map in the physical basis. It attumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, Distribution> PhysMap;
    PhysMap.insert({0, PlusMinusMap.at(0)});


    // Fill in map in the physical basis.
    for (int i = 1; i <= 6; i++)
      {
        PhysMap.insert({i,  ( PlusMinusMap.at(2*i-1) + PlusMinusMap.at(2*i) ) / 2});
        PhysMap.insert({-i, ( PlusMinusMap.at(2*i-1) - PlusMinusMap.at(2*i) ) / 2});
      }

    return PhysMap;
  }
}
