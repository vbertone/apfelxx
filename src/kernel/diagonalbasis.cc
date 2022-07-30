//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#include "apfel/diagonalbasis.h"

namespace apfel
{
  //_________________________________________________________________________________
  DiagonalBasis::DiagonalBasis(int const& nf, int const& offset):
    ConvolutionMap{"DiagonalBasis_" + (std::string) (offset >= 0 ? "p" : "m") + std::to_string(std::abs(offset)) + "_" + std::to_string(nf)}
  {
    for (int k = offset; k < offset + nf; k++)
      _rules[k] = { {k, k, 1} };
  }
}
