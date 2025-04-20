//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/evolutionbasisqcdqed.h"

namespace apfel
{
  //_________________________________________________________________________
  EvolutionBasisQCDQED::EvolutionBasisQCDQED(int const& nd, int const& nu, int const& nl):
    ConvolutionMap{"EvolutionBasisQCDQED_nd" + std::to_string(nd) + "_nu" + std::to_string(nu) + "_nl" + std::to_string(nl)}
  {
    // Helper vectors
    std::vector<int> ActiveDownPlus(nd);
    std::vector<int> ActiveDownMinus(nd);
    std::vector<int> ActiveUpPlus(nu);
    std::vector<int> ActiveUpMinus(nu);
    std::vector<int> ActiveLeptPlus(nl);
    std::vector<int> ActiveLeptMinus(nl);
    for (int id = 0; id < nd; id++)
      {
        ActiveDownPlus[id]  = 8  - id;
        ActiveDownMinus[id] = 11 + id;
      }
    for (int iu = 0; iu < nu; iu++)
      {
        ActiveUpPlus[iu]  = 5  - iu;
        ActiveUpMinus[iu] = 14 + iu;
      }
    for (int il = 0; il < nl; il++)
      {
        ActiveLeptPlus[il]  = 2  - il;
        ActiveLeptMinus[il] = 17 + il;
      }

    // Inactive lepton plus distributions
    for (int i = 2 - nl; i >= 0; i--)
      _rules[i] = {{PPV, i, 0}};

    // Inactive up-type plus distributions
    for (int i = 5 - nu; i > 2; i--)
      _rules[i] = {{PPV, i, 0}};

    // Inactive down-type plus distributions
    for (int i = 8 - nd; i > 5; i--)
      _rules[i] = {{PPV, i, 0}};

    // Inactive down-type minus distributions
    for (int i = 11 + nd; i < 14; i++)
      _rules[i] = {{PPV, i, 0}};

    // Inactive up-type minus distributions
    for (int i = 14 + nu; i < 17; i++)
      _rules[i] = {{PPV, i, 0}};

    // Inactive lepton minus distributions
    for (int i = 17 + nl; i < 20; i++)
      _rules[i] = {{PPV, i, 0}};

    // Active lepton plus distributions
    for (int i : ActiveLeptPlus)
      {
        _rules[i] = {{PPLL, i, 1}};
        for (int j : ActiveDownPlus)
          _rules[i].push_back({PPSLD, j, 1});
        for (int j : ActiveUpPlus)
          _rules[i].push_back({PPSLU, j, 1});
        for (int j : ActiveLeptPlus)
          _rules[i].push_back({PPSLL, j, 1});
        _rules[i].push_back({PLgm, PHOTON, 1});
      }

    // Active up-type plus distributions
    for (int i : ActiveUpPlus)
      {
        _rules[i] = {{PPUU, i, 1}};
        for (int j : ActiveDownPlus)
          _rules[i].push_back({PPSUD, j, 1});
        for (int j : ActiveUpPlus)
          _rules[i].push_back({PPSUU, j, 1});
        for (int j : ActiveLeptPlus)
          _rules[i].push_back({PPSUL, j, 1});
        _rules[i].push_back({PUg, GLUON, 1});
        _rules[i].push_back({PUgm, PHOTON, 1});
      }

    // Active down-type plus distributions
    for (int i : ActiveDownPlus)
      {
        _rules[i] = {{PPDD, i, 1}};
        for (int j : ActiveDownPlus)
          _rules[i].push_back({PPSDD, j, 1});
        for (int j : ActiveUpPlus)
          _rules[i].push_back({PPSDU, j, 1});
        for (int j : ActiveLeptPlus)
          _rules[i].push_back({PPSDL, j, 1});
        _rules[i].push_back({PDg, GLUON, 1});
        _rules[i].push_back({PDgm, PHOTON, 1});
      }

    // Gluon
    _rules[GLUON] = {{Pgg, GLUON, 1}};
    _rules[GLUON].push_back({Pggm, PHOTON, 1});
    for (int i : ActiveDownPlus)
      _rules[GLUON].push_back({PgD, i, 1});
    for (int i : ActiveUpPlus)
      _rules[GLUON].push_back({PgU, i, 1});

    // Photon
    _rules[PHOTON] = {{Pgmgm, PHOTON, 1}};
    _rules[PHOTON].push_back({Pgmg, GLUON, 1});
    for (int i : ActiveDownPlus)
      _rules[PHOTON].push_back({PgmD, i, 1});
    for (int i : ActiveUpPlus)
      _rules[PHOTON].push_back({PgmU, i, 1});
    for (int i : ActiveLeptPlus)
      _rules[PHOTON].push_back({PgmL, i, 1});

    // Active down-type minus distributions
    for (int i : ActiveDownMinus)
      {
        _rules[i] = {{PMDD, i, 1}};
        for (int j : ActiveDownMinus)
          _rules[i].push_back({PPV, j, 1});
        for (int j : ActiveUpMinus)
          _rules[i].push_back({PPV, j, 1});
      }

    // Active up-type minus distributions
    for (int i : ActiveUpMinus)
      {
        _rules[i] = {{PMUU, i, 1}};
        for (int j : ActiveDownMinus)
          _rules[i].push_back({PPV, j, 1});
        for (int j : ActiveUpMinus)
          _rules[i].push_back({PPV, j, 1});
      }

    // Active leptom minus distributions
    for (int i : ActiveLeptMinus)
      _rules[i] = {{PMLL, i, 1}};
  }
}
