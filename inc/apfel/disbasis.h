//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"
#include "apfel/tools.h"

#include <vector>
#include <numeric>

using namespace std;

namespace apfel
{
  /**
   * @brief ConvolutionMap derivation for the QCD evolution basis.
   *
   * This class, following the derivation procedure from
   * ConvolutionMap implements the Basis enumerator with custom tags
   * for the objects.
   */
  class DISNCBasis: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enums
     */
    enum Operand: int {CNS, CS, CG};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @brief The class constructor fot the k-th structure function with no nf dependence
     */
  DISNCBasis(int const& k):
    ConvolutionMap{"DISNCBasis_" + std::to_string(k)}
    {
      // Gluon
      _rules[GLUON] = { {CG, GLUON, 1} };
      // Singlet
      _rules[SIGMA] = { {CS, SIGMA, 1./6.} };
      // Total Valence
      _rules[VALENCE] = { {CS, VALENCE, 1./6.} };
      // Non-singlet distributions
      for (int i = 2; i <= 6; i++)
	{
	  double coef = 0;
	  if (i == k)
	    coef = - 1. / i;
	  else if (i >= k+1)
	    coef = 1. / i / ( i - 1 );

	  // Change sign to T3 and V3
	  if (i == 2)
	    coef *= - 1;

	  _rules[2*i-1] = { {CNS, 2*i-1, coef} };
	  _rules[2*i]   = { {CNS, 2*i,   coef} };
	}
    };

    /**
     * @brief The class constructor fot the total structure function independent of "nf"
     */
  DISNCBasis(vector<double> const& Ch):
    ConvolutionMap{"DISNCBasis_tot"}
    {
      if (Ch.size() != 6)
	throw runtime_exception("DISNCBasis", "The charge vector must have 6 entries.");

      // Sum of the fist nf charges
      const double SumCh = accumulate(Ch.begin(), Ch.end(), 0.);

      // Gluon
      _rules[GLUON] = { {CG, GLUON, SumCh} };
      // Singlet
      _rules[SIGMA] = { {CS, SIGMA, SumCh/6} };
      // Total Valence
      _rules[VALENCE] = { {CS, VALENCE, SumCh/6} };
      // Non-singlet distributions
      for (int j = 2; j <= 6; j++)
	{
	  double coef = 0;
	  for (int i = 1; i <= j; i++)
	    if (i < j)
	      coef += Ch[i-1];
	    else
	      coef += Ch[i-1] * ( 1 - j );
	  coef /= j * ( j - 1 );

	  // Change sign to T3 and V3
	  if (j == 2)
	    coef *= - 1;

	  _rules[2*j-1] = { {CNS, 2*j-1, coef} };
	  _rules[2*j]   = { {CNS, 2*j,   coef} };
	}
    };
  };

}
