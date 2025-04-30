//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/rotations.h"
#include "apfel/evolutionbasisqcd.h"

namespace apfel
{
  //_____________________________________________________________________________
  std::map<int, double> PhysToQCDEv(std::map<int, double> const& InPhysMap)
  {
    // Call function in the physical basis
    std::map<int, double> PhysMap = InPhysMap;

    // Fill in keys that do not exist. Start with the gluon (assume
    // that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap[0] = PhysMap[21];

    // Quarks (fill in with zero if they do not exist)
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap[i] = 0;

    // Fill in map in the QCD evolution basis. It assumes that the
    // gluon has key zero and all keys from -6 to 6 exist.
    std::map<int, double> QCDEvMap;
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
  std::map<int, Distribution> PhysToQCDEv(std::map<int, Distribution> const& InPhysMap)
  {
    // Zero distribution
    const Distribution ZeroDist{InPhysMap.begin()->second.GetGrid(), [] (double const&) -> double { return 0; }};

    // Call function in the physical basis
    std::map<int, Distribution> PhysMap = InPhysMap;

    // Fill in keys that do not exist. Start with the gluon (assume
    // that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap.insert({0, PhysMap.at(21)});

    // Quarks (fill in with zero if they do not exist)
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap.insert({i, ZeroDist});

    // Fill in map in the QCD evolution basis. It assumes that the
    // gluon has key zero and all keys from -6 to 6 exist.
    std::map<int, Distribution> QCDEvMap;
    QCDEvMap.insert({0, PhysMap.at(0)});
    QCDEvMap.insert({1,
                     PhysMap.at(1) + PhysMap.at(-1)
                     + PhysMap.at(2) + PhysMap.at(-2)
                     + PhysMap.at(3) + PhysMap.at(-3)
                     + PhysMap.at(4) + PhysMap.at(-4)
                     + PhysMap.at(5) + PhysMap.at(-5)
                     + PhysMap.at(6) + PhysMap.at(-6)});
    QCDEvMap.insert({2,
                     PhysMap.at(1) - PhysMap.at(-1)
                     + PhysMap.at(2) - PhysMap.at(-2)
                     + PhysMap.at(3) - PhysMap.at(-3)
                     + PhysMap.at(4) - PhysMap.at(-4)
                     + PhysMap.at(5) - PhysMap.at(-5)
                     + PhysMap.at(6) - PhysMap.at(-6)});
    QCDEvMap.insert({3,
                     PhysMap.at(2) + PhysMap.at(-2)
                     - ( PhysMap.at(1) + PhysMap.at(-1) )});
    QCDEvMap.insert({4,
                     PhysMap.at(2) - PhysMap.at(-2)
                     - ( PhysMap.at(1) - PhysMap.at(-1) )});
    QCDEvMap.insert({5,
                     PhysMap.at(1) + PhysMap.at(-1)
                     + PhysMap.at(2) + PhysMap.at(-2)
                     - 2 * ( PhysMap.at(3) + PhysMap.at(-3) )});
    QCDEvMap.insert({6,
                     PhysMap.at(1) - PhysMap.at(-1)
                     + PhysMap.at(2) - PhysMap.at(-2)
                     - 2 * ( PhysMap.at(3) - PhysMap.at(-3) )});
    QCDEvMap.insert({7,
                     PhysMap.at(1) + PhysMap.at(-1)
                     + PhysMap.at(2) + PhysMap.at(-2)
                     + PhysMap.at(3) + PhysMap.at(-3)
                     - 3 * ( PhysMap.at(4) + PhysMap.at(-4) )});
    QCDEvMap.insert({8,
                     PhysMap.at(1) - PhysMap.at(-1)
                     + PhysMap.at(2) - PhysMap.at(-2)
                     + PhysMap.at(3) - PhysMap.at(-3)
                     - 3 * ( PhysMap.at(4) - PhysMap.at(-4) )});
    QCDEvMap.insert({9,
                     PhysMap.at(1) + PhysMap.at(-1)
                     + PhysMap.at(2) + PhysMap.at(-2)
                     + PhysMap.at(3) + PhysMap.at(-3)
                     + PhysMap.at(4) + PhysMap.at(-4)
                     - 4 * ( PhysMap.at(5) + PhysMap.at(-5) )});
    QCDEvMap.insert({10,
                     PhysMap.at(1) - PhysMap.at(-1)
                     + PhysMap.at(2) - PhysMap.at(-2)
                     + PhysMap.at(3) - PhysMap.at(-3)
                     + PhysMap.at(4) - PhysMap.at(-4)
                     - 4 * ( PhysMap.at(5) - PhysMap.at(-5) )});
    QCDEvMap.insert({11,
                     PhysMap.at(1) + PhysMap.at(-1)
                     + PhysMap.at(2) + PhysMap.at(-2)
                     + PhysMap.at(3) + PhysMap.at(-3)
                     + PhysMap.at(4) + PhysMap.at(-4)
                     + PhysMap.at(5) + PhysMap.at(-5)
                     - 5 * ( PhysMap.at(6) + PhysMap.at(-6) )});
    QCDEvMap.insert({12,
                     PhysMap.at(1) - PhysMap.at(-1)
                     + PhysMap.at(2) - PhysMap.at(-2)
                     + PhysMap.at(3) - PhysMap.at(-3)
                     + PhysMap.at(4) - PhysMap.at(-4)
                     + PhysMap.at(5) - PhysMap.at(-5)
                     - 5 * ( PhysMap.at(6) - PhysMap.at(-6) )});

    return QCDEvMap;
  }

  //_____________________________________________________________________________
  Set<Distribution> PhysToQCDEv(std::map<int, Distribution> const& InPhysMap, int const& nf)
  {
    return Set<Distribution> {EvolutionBasisQCD{nf}, PhysToQCDEv(InPhysMap)};
  }

  //_____________________________________________________________________________
  Set<Operator> PhysToQCDEv(std::map<int, Operator> const& InPhysMap, int const& nf)
  {
    // Fill in map in the QCD evolution basis. It assumes that the
    // gluon has key zero and all keys from -6 to 6 exist.
    std::map<int, Operator> QCDEvMap;
    QCDEvMap.insert({0, InPhysMap.at(0)});
    QCDEvMap.insert({1,
                     InPhysMap.at(1) + InPhysMap.at(-1)
                     + InPhysMap.at(2) + InPhysMap.at(-2)
                     + InPhysMap.at(3) + InPhysMap.at(-3)
                     + InPhysMap.at(4) + InPhysMap.at(-4)
                     + InPhysMap.at(5) + InPhysMap.at(-5)
                     + InPhysMap.at(6) + InPhysMap.at(-6)});
    QCDEvMap.insert({2,
                     InPhysMap.at(1) - InPhysMap.at(-1)
                     + InPhysMap.at(2) - InPhysMap.at(-2)
                     + InPhysMap.at(3) - InPhysMap.at(-3)
                     + InPhysMap.at(4) - InPhysMap.at(-4)
                     + InPhysMap.at(5) - InPhysMap.at(-5)
                     + InPhysMap.at(6) - InPhysMap.at(-6)});
    QCDEvMap.insert({3,
                     InPhysMap.at(2) + InPhysMap.at(-2)
                     - ( InPhysMap.at(1) + InPhysMap.at(-1) )});
    QCDEvMap.insert({4,
                     InPhysMap.at(2) - InPhysMap.at(-2)
                     - ( InPhysMap.at(1) - InPhysMap.at(-1) )});
    QCDEvMap.insert({5,
                     InPhysMap.at(1) + InPhysMap.at(-1)
                     + InPhysMap.at(2) + InPhysMap.at(-2)
                     - 2 * ( InPhysMap.at(3) + InPhysMap.at(-3) )});
    QCDEvMap.insert({6,
                     InPhysMap.at(1) - InPhysMap.at(-1)
                     + InPhysMap.at(2) - InPhysMap.at(-2)
                     - 2 * ( InPhysMap.at(3) - InPhysMap.at(-3) )});
    QCDEvMap.insert({7,
                     InPhysMap.at(1) + InPhysMap.at(-1)
                     + InPhysMap.at(2) + InPhysMap.at(-2)
                     + InPhysMap.at(3) + InPhysMap.at(-3)
                     - 3 * ( InPhysMap.at(4) + InPhysMap.at(-4) )});
    QCDEvMap.insert({8,
                     InPhysMap.at(1) - InPhysMap.at(-1)
                     + InPhysMap.at(2) - InPhysMap.at(-2)
                     + InPhysMap.at(3) - InPhysMap.at(-3)
                     - 3 * ( InPhysMap.at(4) - InPhysMap.at(-4) )});
    QCDEvMap.insert({9,
                     InPhysMap.at(1) + InPhysMap.at(-1)
                     + InPhysMap.at(2) + InPhysMap.at(-2)
                     + InPhysMap.at(3) + InPhysMap.at(-3)
                     + InPhysMap.at(4) + InPhysMap.at(-4)
                     - 4 * ( InPhysMap.at(5) + InPhysMap.at(-5) )});
    QCDEvMap.insert({10,
                     InPhysMap.at(1) - InPhysMap.at(-1)
                     + InPhysMap.at(2) - InPhysMap.at(-2)
                     + InPhysMap.at(3) - InPhysMap.at(-3)
                     + InPhysMap.at(4) - InPhysMap.at(-4)
                     - 4 * ( InPhysMap.at(5) - InPhysMap.at(-5) )});
    QCDEvMap.insert({11,
                     InPhysMap.at(1) + InPhysMap.at(-1)
                     + InPhysMap.at(2) + InPhysMap.at(-2)
                     + InPhysMap.at(3) + InPhysMap.at(-3)
                     + InPhysMap.at(4) + InPhysMap.at(-4)
                     + InPhysMap.at(5) + InPhysMap.at(-5)
                     - 5 * ( InPhysMap.at(6) + InPhysMap.at(-6) )});
    QCDEvMap.insert({12,
                     InPhysMap.at(1) - InPhysMap.at(-1)
                     + InPhysMap.at(2) - InPhysMap.at(-2)
                     + InPhysMap.at(3) - InPhysMap.at(-3)
                     + InPhysMap.at(4) - InPhysMap.at(-4)
                     + InPhysMap.at(5) - InPhysMap.at(-5)
                     - 5 * ( InPhysMap.at(6) - InPhysMap.at(-6) )});

    return Set<Operator> {EvolutionOperatorBasisQCD{nf}, QCDEvMap};
  }

  //_____________________________________________________________________________
  std::map<int, double> QCDEvToPhys(std::map<int, double> const& QCDEvMap)
  {
    // Fill in map in the physical basis. It assumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, double> PhysMap;
    PhysMap[0]  = QCDEvMap.at(0);
    PhysMap[21] = QCDEvMap.at(0);

    // Perform the rotation
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
    // Fill in map in the physical basis. It assumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, Distribution> PhysMap;
    PhysMap.insert({0, QCDEvMap.at(0)});

    // Perform the rotation
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
    // Fill in map in the physical basis. It assumes that the gluon
    // has key zero and all keys from 0 to 12 exist.
    std::map<int, Operator> PhysMap;
    PhysMap.insert({0, QCDEvMap.at(0)});

    // Perform the rotation
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
    // Call function in the physical basis
    std::map<int, double> PhysMap = InPhysMap;

    // Fill in keys that do not exist. Start with the gluon (assumes
    // that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap[0] = PhysMap[21];

    // Quarks (fill in with zero if they do not exist).
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap[i] = 0;

    // Fill in map in the PlusMinus basis. It assumes that the gluon
    // has key zero and all keys from -6 to 6 exist.
    std::map<int, double> PlusMinusMap;
    PlusMinusMap[0 + 6] = PhysMap.at(0);
    for (int i = 1; i <= 6; i++)
      {
        PlusMinusMap[  i + 6] = PhysMap.at(i) + PhysMap.at(-i);
        PlusMinusMap[- i + 6] = PhysMap.at(i) - PhysMap.at(-i);
      }
    return PlusMinusMap;
  }

  //_____________________________________________________________________________
  std::map<int, Distribution> PhysToPlusMinus(std::map<int, Distribution> const& InPhysMap)
  {
    // Zero distribution
    const Distribution ZeroDist{InPhysMap.begin()->second.GetGrid(), [] (double const&) -> double { return 0; }};

    // Call function in the physical basis
    std::map<int, Distribution> PhysMap = InPhysMap;

    // Fill in keys that do not exist. Start with the gluon (assumes
    // that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap.insert({0, PhysMap.at(21)});

    // Quarks (fill in with zero if they do not exist).
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap.insert({i, ZeroDist});

    // Fill in map in the PlusMinus basis. It assumes that the gluon
    // has key zero and all keys from -6 to 6 exist.
    std::map<int, Distribution> PlusMinusMap;
    PlusMinusMap.insert({0 + 6, PhysMap.at(0)});
    for (int i = 1; i <= 6; i++)
      {
        PlusMinusMap.insert({  i + 6, PhysMap.at(i) + PhysMap.at(-i)});
        PlusMinusMap.insert({- i + 6, PhysMap.at(i) - PhysMap.at(-i)});
      }
    return PlusMinusMap;
  }

  //_____________________________________________________________________________
  std::map<int, double> PlusMinusToPhys(std::map<int, double> const& PlusMinusMap)
  {
    // Fill in map in the physical basis
    std::map<int, double> PhysMap;
    PhysMap[0]  = PlusMinusMap.at(0 + 6);
    PhysMap[21] = PlusMinusMap.at(0 + 6);
    for (int i = 1; i <= 6; i++)
      {
        PhysMap[i]  = ( PlusMinusMap.at(i + 6) + PlusMinusMap.at(- i + 6) ) / 2;
        PhysMap[-i] = ( PlusMinusMap.at(i + 6) - PlusMinusMap.at(- i + 6) ) / 2;
      }
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, Distribution> PlusMinusToPhys(std::map<int, Distribution> const& PlusMinusMap)
  {
    // Fill in map in the physical basis
    std::map<int, Distribution> PhysMap;
    PhysMap.insert({0,  PlusMinusMap.at(0 + 6)});
    for (int i = 1; i <= 6; i++)
      {
        PhysMap.insert({i,  ( PlusMinusMap.at(i + 6) + PlusMinusMap.at(- i + 6) ) / 2});
        PhysMap.insert({-i, ( PlusMinusMap.at(i + 6) - PlusMinusMap.at(- i + 6) ) / 2});
      }
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, double> PlusMinusQCDQEDToPhys(std::map<int, double> const& PlusMinusMap)
  {
    // Fill in map in the physical basis
    std::map<int, double> PhysMap;
    // Gluon
    PhysMap[0]   = PlusMinusMap.at(9);
    PhysMap[21]  = PlusMinusMap.at(9);
    // Photon
    PhysMap[22]  = PlusMinusMap.at(10);
    // Leptons
    PhysMap[15]  = ( PlusMinusMap.at(0) + PlusMinusMap.at(19) ) / 2;
    PhysMap[-15] = ( PlusMinusMap.at(0) - PlusMinusMap.at(19) ) / 2;
    PhysMap[13]  = ( PlusMinusMap.at(1) + PlusMinusMap.at(18) ) / 2;
    PhysMap[-13] = ( PlusMinusMap.at(1) - PlusMinusMap.at(18) ) / 2;
    PhysMap[11]  = ( PlusMinusMap.at(2) + PlusMinusMap.at(17) ) / 2;
    PhysMap[-11] = ( PlusMinusMap.at(2) - PlusMinusMap.at(17) ) / 2;
    // Up-type quarks
    PhysMap[6]   = ( PlusMinusMap.at(3) + PlusMinusMap.at(16) ) / 2;
    PhysMap[-6]  = ( PlusMinusMap.at(3) - PlusMinusMap.at(16) ) / 2;
    PhysMap[4]   = ( PlusMinusMap.at(4) + PlusMinusMap.at(15) ) / 2;
    PhysMap[-4]  = ( PlusMinusMap.at(4) - PlusMinusMap.at(15) ) / 2;
    PhysMap[2]   = ( PlusMinusMap.at(5) + PlusMinusMap.at(14) ) / 2;
    PhysMap[-2]  = ( PlusMinusMap.at(5) - PlusMinusMap.at(14) ) / 2;
    // Down-type quarks
    PhysMap[5]   = ( PlusMinusMap.at(6) + PlusMinusMap.at(13) ) / 2;
    PhysMap[-5]  = ( PlusMinusMap.at(6) - PlusMinusMap.at(13) ) / 2;
    PhysMap[3]   = ( PlusMinusMap.at(7) + PlusMinusMap.at(12) ) / 2;
    PhysMap[-3]  = ( PlusMinusMap.at(7) - PlusMinusMap.at(12) ) / 2;
    PhysMap[1]   = ( PlusMinusMap.at(8) + PlusMinusMap.at(11) ) / 2;
    PhysMap[-1]  = ( PlusMinusMap.at(8) - PlusMinusMap.at(11) ) / 2;
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, Distribution> PlusMinusQCDQEDToPhys(std::map<int, Distribution> const& PlusMinusMap)
  {
    // Fill in map in the physical basis
    std::map<int, Distribution> PhysMap;
    // Gluon
    PhysMap.insert({0,   PlusMinusMap.at(9)});
    // Photon
    PhysMap.insert({22,  PlusMinusMap.at(10)});
    // Leptons
    PhysMap.insert({15,  ( PlusMinusMap.at(0) + PlusMinusMap.at(19) ) / 2});
    PhysMap.insert({-15, ( PlusMinusMap.at(0) - PlusMinusMap.at(19) ) / 2});
    PhysMap.insert({13,  ( PlusMinusMap.at(1) + PlusMinusMap.at(18) ) / 2});
    PhysMap.insert({-13, ( PlusMinusMap.at(1) - PlusMinusMap.at(18) ) / 2});
    PhysMap.insert({11,  ( PlusMinusMap.at(2) + PlusMinusMap.at(17) ) / 2});
    PhysMap.insert({-11, ( PlusMinusMap.at(2) - PlusMinusMap.at(17) ) / 2});
    // Up-type quarks
    PhysMap.insert({6,   ( PlusMinusMap.at(3) + PlusMinusMap.at(16) ) / 2});
    PhysMap.insert({-6,  ( PlusMinusMap.at(3) - PlusMinusMap.at(16) ) / 2});
    PhysMap.insert({4,   ( PlusMinusMap.at(4) + PlusMinusMap.at(15) ) / 2});
    PhysMap.insert({-4,  ( PlusMinusMap.at(4) - PlusMinusMap.at(15) ) / 2});
    PhysMap.insert({2,   ( PlusMinusMap.at(5) + PlusMinusMap.at(14) ) / 2});
    PhysMap.insert({-2,  ( PlusMinusMap.at(5) - PlusMinusMap.at(14) ) / 2});
    // Down-type quarks
    PhysMap.insert({5,   ( PlusMinusMap.at(6) + PlusMinusMap.at(13) ) / 2});
    PhysMap.insert({-5,  ( PlusMinusMap.at(6) - PlusMinusMap.at(13) ) / 2});
    PhysMap.insert({3,   ( PlusMinusMap.at(7) + PlusMinusMap.at(12) ) / 2});
    PhysMap.insert({-3,  ( PlusMinusMap.at(7) - PlusMinusMap.at(12) ) / 2});
    PhysMap.insert({1,   ( PlusMinusMap.at(8) + PlusMinusMap.at(11) ) / 2});
    PhysMap.insert({-1,  ( PlusMinusMap.at(8) - PlusMinusMap.at(11) ) / 2});
    return PhysMap;
  }

  //_____________________________________________________________________________
  std::map<int, double> PhysToPlusMinusQCDQED(std::map<int, double> const& InPhysMap)
  {
    // Call function in the physical basis
    std::map<int, double> PhysMap = InPhysMap;

    // Fill in keys that do not exist. Start with the gluon (assumes
    // that the ID is 21).
    if (PhysMap.find(0) == PhysMap.end())
      PhysMap[0] = PhysMap[21];

    // Photon (fill in with zero if they do not exist).
    if (PhysMap.find(22) == PhysMap.end())
      PhysMap[22] = 0;

    // Quarks (fill in with zero if they do not exist).
    for (int i = -6; i <= 6; i++)
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap[i] = 0;

    // Leptons (fill in with zero if they do not exist).
    for (int i : std::vector<int> {-15, -13, -11, 11, 13, 15})
      if (PhysMap.find(i) == PhysMap.end())
        PhysMap[i] = 0;

    // Fill in map in the PlusMinusQCDQED basis
    std::map<int, double> PlusMinusMapQCDQED;
    PlusMinusMapQCDQED[0]  = PhysMap.at(15) + PhysMap.at(-15);
    PlusMinusMapQCDQED[1]  = PhysMap.at(13) + PhysMap.at(-13);
    PlusMinusMapQCDQED[2]  = PhysMap.at(11) + PhysMap.at(-11);
    PlusMinusMapQCDQED[3]  = PhysMap.at(6) + PhysMap.at(-6);
    PlusMinusMapQCDQED[4]  = PhysMap.at(4) + PhysMap.at(-4);
    PlusMinusMapQCDQED[5]  = PhysMap.at(2) + PhysMap.at(-2);
    PlusMinusMapQCDQED[6]  = PhysMap.at(5) + PhysMap.at(-5);
    PlusMinusMapQCDQED[7]  = PhysMap.at(3) + PhysMap.at(-3);
    PlusMinusMapQCDQED[8]  = PhysMap.at(1) + PhysMap.at(-1);
    PlusMinusMapQCDQED[9]  = PhysMap.at(0);
    PlusMinusMapQCDQED[10] = PhysMap.at(22);
    PlusMinusMapQCDQED[11] = PhysMap.at(1) - PhysMap.at(-1);
    PlusMinusMapQCDQED[12] = PhysMap.at(3) - PhysMap.at(-3);
    PlusMinusMapQCDQED[13] = PhysMap.at(5) - PhysMap.at(-5);
    PlusMinusMapQCDQED[14] = PhysMap.at(2) - PhysMap.at(-2);
    PlusMinusMapQCDQED[15] = PhysMap.at(4) - PhysMap.at(-4);
    PlusMinusMapQCDQED[16] = PhysMap.at(6) - PhysMap.at(-6);
    PlusMinusMapQCDQED[17] = PhysMap.at(11) - PhysMap.at(-11);
    PlusMinusMapQCDQED[18] = PhysMap.at(13) - PhysMap.at(-13);
    PlusMinusMapQCDQED[19] = PhysMap.at(15) - PhysMap.at(-15);
    return PlusMinusMapQCDQED;
  }
}
