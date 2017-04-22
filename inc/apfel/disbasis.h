//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

#include <vector>
#include <numeric>

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
     * @brief The class constructor fot the k-th structure function
     */
  DISNCBasis(int const& k, int const& nf):
    ConvolutionMap{"DISNCBasis_" + std::to_string(k) + "_" + std::to_string(nf)}
    {
      // Gluon
      _rules[GLUON] = { {CG, GLUON, 1} };
      // Singlet/total valence
      _rules[DTOT] = { {CT, DTOT, 1./nf} };
      // Non-singlet distributions
      for (int i = D3; i <= D35; i++)
	{
	  double coef = 0;
	  if (i == k)
	    coef = - 1. / i;
	  else if (i >= k+1 && i <= nf)
	    coef = 1. / i / ( i - 1 );
	  _rules[i] = { {CNS, i, coef} };
	}

      // Set all multiplicative coefficients to zero if "k" > "nf".
      // This is equivalent to theta(Q^2 - mh^2).
      if (k > nf)
	for (int i = GLUON; i <= D35; i++)
	  _rules[i] = { {CNS, i, 0} };
    };


    /**
     * @brief The class constructor fot the total structure function
     */
  DISNCBasis(vector<double> const& Ch, int const& nf):
    ConvolutionMap{"DISNCBasis_tot_" + std::to_string(nf)}
    {
      if (Ch.size() != 6)
	throw runtime_exception("DISNCBasis", "The charge vector must have 6 entries.");

      // Sum of the fist nf charges
      const double SumCh = accumulate(Ch.begin(), Ch.begin() + nf, 0.);

      // Gluon
      _rules[GLUON] = { {CG, GLUON, SumCh} };
      // Singlet/total valence
      _rules[DTOT] = { {CT, DTOT, SumCh/nf} };
      // Non-singlet distributions
      for (int j = 2; j <= 6; j++)
	{
	  double coef = 0;
	  if (j <= nf)
	    {
	      for (int i = 1; i <= j; i++)
		if (i < j)
		  coef += Ch[i-1];
		else
		  coef += Ch[i-1] * ( 1 - j );
	      coef /= j * ( j - 1 );
	    }
	  _rules[j] = { {CNS, j, coef} };
	}
    };
  };

}
