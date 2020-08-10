//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/evolutionbasisqcd.h"

namespace apfel
{
  //_________________________________________________________________________
  EvolutionBasisQCD::EvolutionBasisQCD(int const& nf):
    ConvolutionMap{"EvolutionBasisQCD_" + std::to_string(nf)}
  {
    _rules[GLUON] = { {PGG, GLUON, 1}, {PGQ, SIGMA, 1} };
    for (int j = nf + 1; j <= 6; j++)
      _rules[GLUON].push_back({PGQ, 2 * j - 1, 6. / j / ( j - 1 )});

    _rules[SIGMA]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
    for (int j = nf + 1; j <= 6; j++)
      _rules[SIGMA].push_back({PQQ, 2 * j - 1, 6. / j / ( j - 1 )});

    _rules[VALENCE] = { {PNSV, VALENCE, 1} };

    for (int i = 2; i <= nf; i++)
      {
        _rules[2 * i - 1] = { {PNSP, 2 * i - 1, 1} };
        _rules[2 * i]     = { {PNSM, 2 * i,     1} };
      }

    for (int i = nf + 1; i <= 6; i++)
      {
        _rules[2 * i - 1]   = { {PQG, GLUON, 1}, {PQQ, SIGMA, 1} };
        for (int j = nf + 1; j <= 6; j++)
          _rules[2 * i - 1].push_back({PQQ, 2 * j - 1, 6. / j / ( j - 1 )});
        _rules[2 * i]     = { {PNSV, 2 * i, 1} };
      }
  }

  //_________________________________________________________________________
  EvolutionOperatorBasisQCD::EvolutionOperatorBasisQCD(int const& nf):
    ConvolutionMap{"EvolutionOperatorBasisQCD_" + std::to_string(nf)}
  {
    // Allocate EvolutionBasisQCD object to retrieve the splitting
    // matrix rules
    const EvolutionBasisQCD eb{nf};

    // Get matrix of coefficients
    const matrix<std::vector<double>> rc = eb.GetRuleMatrix();

    // Get matrix of operator indices
    const matrix<std::vector<int>> ri = eb.GetRuleIndices();

    // Now construct set of rules
    for (int i = 0; i < 13; i++)
      for (int j = 0; j < 13; j++)
        for (int k = 0; k < 13; k++)
          {
            if (rc(i, k).empty() || Gkj.count({k, j}) == 0)
              continue;

            for (int l = 0; l < (int) rc(i, k).size(); l++)
              _rules[Gkj.at({i, j})].push_back({ri(i, k)[l], Gkj.at({k, j}), rc(i, k)[l]});
          }
  }

  //_________________________________________________________________________
  EvolveDistributionsBasisQCD::EvolveDistributionsBasisQCD():
    ConvolutionMap{"EvolveDistributionsBasisQCD"}
  {
    // Construct set of rules
    for (int k = 0; k < 13; k++)
      for (int j = 0; j < 13; j++)
        {
          if (Gkj.count({k, j}) == 0)
            continue;
          _rules[k].push_back({Gkj.at({k, j}), j, 1});
        }
  }
}
