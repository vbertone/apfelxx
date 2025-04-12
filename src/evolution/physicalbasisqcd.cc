//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/physicalbasisqcd.h"

namespace apfel
{
  //_________________________________________________________________________
  PhysicalBasisQCD::PhysicalBasisQCD(int const& nf):
    ConvolutionMap{"PhysicalBasisQCD_" + std::to_string(nf)}
  {
    // Minus-type distributions
    for (int i = - 6 + 6; i <= - nf - 1 + 6; i++)
      _rules[i] = {{PPV, i, 0}};
    for (int i = - nf + 6; i <= - 1 + 6; i++)
      for (int j = - nf + 6; j <= - 1 + 6; j++)
        if (j == i)
          _rules[i].push_back({PNV, j, 1});
        else
          _rules[i].push_back({PPV, j, 1});

    // Gluon
    _rules[GLUON] = {{PGG, GLUON, 1}};
    for (int j = 1 + 6; j <= nf + 6; j++)
      _rules[GLUON].push_back({PGQ, j, 1});

    // Plus-type distributions
    for (int i = 1 + 6; i <= nf + 6; i++)
      {
        _rules[i] = {{PQG, 0 + 6, 1}};
        for (int j = 1 + 6; j <= nf + 6; j++)
          if (j == i)
            _rules[i].push_back({PNS, j, 1});
          else
            _rules[i].push_back({PPS, j, 1});
      }
    for (int i = nf + 1 + 6; i <= 6 + 6; i++)
      _rules[i] = {{PPS, i, 0}};
  }

  //_________________________________________________________________________
  PhysicalOperatorBasisQCD::PhysicalOperatorBasisQCD(int const& nf):
    ConvolutionMap{"PhysicalOperatorBasisQCD_" + std::to_string(nf)}
  {
    // Allocate PhysicalBasisQCD object to retrieve the splitting
    // matrix rules.
    const PhysicalBasisQCD eb{nf};

    // Get matrix of coefficients
    const matrix<std::vector<double>> rc = eb.GetRuleMatrix();

    // Get matrix of operator indices
    const matrix<std::vector<int>> ri = eb.GetRuleIndices();

    // Now construct set of rules
    for (int i = 0; i < 13; i++)
      for (int j = 0; j < 13; j++)
        for (int k = 0; k < 13; k++)
          {
            if (rc(i, k).empty() || GkjPhys.count({k, j}) == 0)
              continue;

            for (int l = 0; l < (int) rc(i, k).size(); l++)
              _rules[GkjPhys.at({i, j})].push_back({ri(i, k)[l], GkjPhys.at({k, j}), rc(i, k)[l]});
          }
  }

  //_________________________________________________________________________
  PhysicalEvolveDistributionsBasisQCD::PhysicalEvolveDistributionsBasisQCD():
    ConvolutionMap{"PhysicalEvolveDistributionsBasisQCD"}
  {
    // Construct set of rules
    for (int k = 0; k < 13; k++)
      for (int j = 0; j < 13; j++)
        {
          if (GkjPhys.count({k, j}) == 0)
            continue;
          _rules[k].push_back({GkjPhys.at({k, j}), j, 1});
        }
  }
}
