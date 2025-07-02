//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/matchingbasisqcdqed.h"
#include "apfel/evolutionbasisqcdqed.h"

namespace apfel
{
  //_________________________________________________________________________________
  MatchingBasisQCDQED::MatchingBasisQCDQED(int const& nd, int const& nu, int const& nl, PartonSpecies const& species):
    ConvolutionMap{"MatchingBasisQCDQED_nd" + std::to_string(nd) + "_nu" + std::to_string(nu) + "_nl" + std::to_string(nl) + "_species"  + std::to_string(species)}
  {
    // Leading-order matching
    for (int i = Object::TAUP; i <= Object::TAUM; i++)
      _rules[i] = {{ONE, i, 1}};

    // Helper vectors
    std::vector<int> LightQuarksPlus(nd + nu);
    std::vector<int> LightQuarksMinus(nd + nu);
    for (int id = 0; id < nd; id++)
      {
        LightQuarksPlus[id]  = 8  - id;
        LightQuarksMinus[id] = 11 + id;
      }
    for (int iu = 0; iu < nu; iu++)
      {
        LightQuarksPlus[nd + iu]  = 5  - iu;
        LightQuarksMinus[nd + iu] = 14 + iu;
      }

    // Heavy-parton plus-type index
    int ih = -1;
    if (species == PartonSpecies::DOWNTYPEQUARK)
      ih = 8 - nd;
    else if (species == PartonSpecies::UPTYPEQUARK)
      ih = 5 - nu;
    else if (species == PartonSpecies::CHARGEDLEPTON)
      ih = 2 - nl;
    else
      return;

    // Heavy parton up-type
    _rules[ih].push_back({KQg, GLUON, 1});
    for (int j : LightQuarksPlus)
      _rules[ih].push_back({KQqp, j, 1});
    _rules[ih].push_back({KXX, ih, 1});

    // Light-quarks up-type distributions
    for (int i : LightQuarksPlus)
      {
        _rules[i].push_back({Kqg, GLUON, 1});
        _rules[i].push_back({KNSq, i, 1});
        for (int j : LightQuarksPlus)
          _rules[i].push_back({Kqqp, j, 1});
      }

    // Gluon
    _rules[GLUON].push_back({Kgg, GLUON, 1});
    for (int j : LightQuarksPlus)
      _rules[GLUON].push_back({Kgq, j, 1});
    _rules[GLUON].push_back({KgQ, ih, 1});

    // Photon
    _rules[PHOTON].push_back({KXgmgm, PHOTON, 1});
    _rules[PHOTON].push_back({KgmX, ih, 1});

    // Light-quarks minus-type distributions
    for (int i : LightQuarksMinus)
      _rules[i].push_back({KNSqm, i, 1});
  }

  //_________________________________________________________________________________
  MatchingOperatorBasisQCDQED::MatchingOperatorBasisQCDQED(int const& nd, int const& nu, int const& nl, PartonSpecies const& species):
    ConvolutionMap{"MatchingOperatorBasisQCDQED_nd" + std::to_string(nd) + "_nu" + std::to_string(nu) + "_nl" + std::to_string(nl) + "_species"  + std::to_string(species)}
  {
    // Allocate MatchingBasisQCD object to retrieve the splitting
    // matrix rules.
    const MatchingBasisQCDQED mb{nd, nu, nl, species};

    // Get matrix of coefficients
    const matrix<std::vector<double>> rc = mb.GetRuleMatrix();

    // Get matrix of operator indices
    const matrix<std::vector<int>> ri = mb.GetRuleIndices();

    // Now construct set of rules
    for (int i = 0; i < 20; i++)
      for (int j = 0; j < 20; j++)
        for (int k = 0; k < 20; k++)
          {
            if (rc(i, k).empty() || GkjQCDQED.count({k, j}) == 0)
              continue;

            for (int l = 0; l < (int) rc(i, k).size(); l++)
              _rules[GkjQCDQED.at({i, j})].push_back({ri(i, k)[l], GkjQCDQED.at({k, j}), rc(i, k)[l]});
          }
  }
}
