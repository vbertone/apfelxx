//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

using namespace std;

namespace apfel
{
  /**
   * @brief ConvolutionMap derivation for the QCD evolution basis.
   *
   * This class, following the derivation procedure from ConvolutionMap
   * implements the Basis enumerator with custom tags for the objects.
   */
  class DISNCBasis: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enums
     */
    enum Operand: int {CNS, CT, CG};
    enum Object:  int {GLUON, DTOT, D3, D8, D15, D24, D35};

    /**
     * @brief The class constructor
     */
  DISNCBasis(int const& k, int const& nf):
    ConvolutionMap{"DISNCBasis_" + std::to_string(k) + "_" + std::to_string(nf)}
    {
      double coeft = 1;
      if (k > nf)
	 coeft = 0;
      // Gluon
      _rules[GLUON] = { {CG, GLUON, coeft} };
      // Singlet/total valence
      _rules[DTOT] = { {CT, DTOT, coeft/nf} };
      // Non-singlet distributions
      for (int i = D3; i <= D35; i++)
	{
	  double coefns = 0;
	  if (i == k)
	    coefns = - 1. / i;
	  else if (i >= k+1 && i <= nf)
	    coefns = 1. / i / ( i - 1 );
	  if (k > nf)
	    coefns = 0;
	  _rules[i] = { {CNS, i, coefns} };
	}
    };
  };

}
