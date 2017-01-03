//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/distribution.h"
#include "apfel/grid.h"
#include "apfel/subgrid.h"

namespace apfel
{
  //_________________________________________________________________________
  Distribution::Distribution(Grid const& gr):
    LagrangeInterpolator{gr}
  {
  }

  //_________________________________________________________________________
  Distribution::Distribution(const Distribution &obj, vector<vector<double>> const& distsubgrid):
    LagrangeInterpolator{obj._grid}
  {
    _distributionSubGrid = distsubgrid;
  }
}
