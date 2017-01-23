//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#include "apfel/basismap.h"

namespace apfel {

  //_________________________________________________________________________________
  BasisMap::BasisMap(string const& name):
    _name(name)
  {
  }

  //_________________________________________________________________________________
  FlvrBasis::FlvrBasis():
    BasisMap{"FlvrBasis"}
  {
    // g = Pgg*g + Sum Pqq*q
    _rules[GLU] = {
                    {PGG, GLU, +1},
                    {PGQ, U, +1},
                    {PGQ, D, +1},
                    {PGQ, S, +1},
                    {PGQ, UBAR, +1},
                    {PGQ, DBAR, +1},
                    {PGQ, SBAR, +1},
                  };

    // q = Pqq*q + Pgg*g
    for (const auto &v: {U,D,S,UBAR,DBAR,SBAR})
    _rules[v] = {
                  {PQQ, v, +1},
                  {PGG, GLU, +1}
                };
  }

}
