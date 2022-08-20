//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/evolbasisplusminus.h"

namespace apfel
{
  //_________________________________________________________________________________
  EvolBasisPlusMinus::EvolBasisPlusMinus():
    ConvolutionMap{"EvolBasisPlusMinus"}
  {
    _rules[GLUON] = { {PGG, GLUON, 1}, {PGD, DWP, 1}, {PGU, UPP, 1}, {PGS, STP, 1}, {PGC, CHP, 1}, {PGB, BTP, 1}, {PGT, TPP, 1}};
    _rules[DWP]   = { {PDD, DWP,   1}, {PPS, UPP, 1}, {PPS, STP, 1}, {PPS, CHP, 1}, {PPS, BTP, 1}, {PPS, TPP, 1}, {PDG, GLUON, 1}};
    _rules[DWM]   = { {PDD, DWM,   1}, {PMP, DWM, 1}};
    _rules[UPP]   = { {PPS, DWP,   1}, {PUU, UPP, 1}, {PPS, STP, 1}, {PPS, CHP, 1}, {PPS, BTP, 1}, {PPS, TPP, 1}, {PUG, GLUON, 1}};
    _rules[UPM]   = { {PUU, UPM,   1}, {PMP, UPM, 1}};
    _rules[STP]   = { {PPS, DWP,   1}, {PPS, UPP, 1}, {PSS, STP, 1}, {PPS, CHP, 1}, {PPS, BTP, 1}, {PPS, TPP, 1}, {PSG, GLUON, 1}};
    _rules[STM]   = { {PSS, STM,   1}, {PMP, STM, 1}};
    _rules[CHP]   = { {PPS, DWP,   1}, {PPS, UPP, 1}, {PPS, STP, 1}, {PCC, CHP, 1}, {PPS, BTP, 1}, {PPS, TPP, 1}, {PCG, GLUON, 1}};
    _rules[CHM]   = { {PCC, CHM,   1}, {PMP, CHM, 1}};
    _rules[BTP]   = { {PPS, DWP,   1}, {PPS, UPP, 1}, {PPS, STP, 1}, {PPS, CHP, 1}, {PBB, BTP, 1}, {PPS, TPP, 1}, {PBG, GLUON, 1}};
    _rules[BTM]   = { {PBB, BTM,   1}, {PMP, BTM, 1}};
    _rules[TPP]   = { {PPS, DWP,   1}, {PPS, UPP, 1}, {PPS, STP, 1}, {PPS, CHP, 1}, {PPS, BTP, 1}, {PTT, TPP, 1}, {PTG, GLUON, 1}};
    _rules[TPM]   = { {PTT, TPM,   1}, {PMP, TPM, 1}};
  }
}
