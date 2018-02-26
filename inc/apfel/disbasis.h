//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"
#include "apfel/messages.h"

#include <vector>
#include <numeric>

using namespace std;

namespace apfel
{
  /**
   * @brief The map enumerators for the operands and the
   * distributions.
   */
  enum Operand: int {CNS, CS, CG};
  enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

  /**
   * @defgroup DISBases DIS convolution maps
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for NC and CC DIS structure functions.
   */
  ///@{
  /**
   * @brief The DISNCBasis class is a derived of ConvolutionMap
   * specialised for the computation of the NC DIS structure
   * functions.
   */
  class DISNCBasis: public ConvolutionMap
  {
  public:
    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DISNCBasis default constructor for the k-th
     * component of the structure functions.
     * @param k: index that identifies the component of the structure function
     * @param fact: factor that multiplies the whole structure function
     */
  DISNCBasis(int const& k, double const& fact = 1):
    ConvolutionMap{"DISNCBasis_" + std::to_string(k)}
    {
      _rules[GLUON]   = { {CG, GLUON, fact} };
      _rules[SIGMA]   = { {CS, SIGMA, fact / 6} };
      _rules[VALENCE] = { {CS, VALENCE, fact / 6} };
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

	  // Multiply the coefficient by the overall factor.
	  coef *= fact;

	  _rules[2*i-1] = { {CNS, 2*i-1, coef} };
	  _rules[2*i]   = { {CNS, 2*i,   coef} };
	}
    };

    /**
     * @brief The DISNCBasis default constructor for the total structure functions.
     * @param Ch: vector of the effective quark charges
     */
  DISNCBasis(vector<double> const& Ch):
    ConvolutionMap{"DISNCBasis_tot"}
    {
      if (Ch.size() != 6)
	throw runtime_error(error("DISNCBasis", "The charge vector must have 6 entries."));

      // Sum of the charges
      const double SumCh = accumulate(Ch.begin(), Ch.end(), 0.);

      _rules[GLUON]   = { {CG, GLUON, SumCh} };
      _rules[SIGMA]   = { {CS, SIGMA, SumCh/6} };
      _rules[VALENCE] = { {CS, VALENCE, SumCh/6} };
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
    ///@}
  };

  /**
   * @brief The DISCCBasis class is a derived of ConvolutionMap
   * specialised for the computation of the CC DIS structure
   * functions.
   */
  class DISCCBasis: public ConvolutionMap
  {
  public:
    /**
     * @brief Map between one single index and the CKM matrix
     * elements:
     *
     * 1 - Vud2
     * 2 - Vus2
     * 3 - Vub2
     * 4 - Vcd2
     * 5 - Vcs2
     * 6 - Vcb2
     * 7 - Vtd2
     * 8 - Vts2
     * 9 - Vtb2
     */
    map<int,pair<int,int>> Vij = {
      {1,{1,1}}, {2,{1,2}}, {3,{1,3}},
      {4,{2,1}}, {5,{2,2}}, {6,{2,3}},
      {7,{3,1}}, {8,{3,2}}, {9,{3,3}} };

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DISCCBasis default constructor for the (i,k)-th
     * component of the structure functions.
     * @param l: index that identifies the component of the structure function
     * @param Is3: switch to tell the constructure whether the structure function is F3 or not
     * @param fact: factor that multiplies the whole structure function
     */
  DISCCBasis(int const& l, bool const& Is3, double const& fact = 1):
    ConvolutionMap{"DISCCBasis_" + to_string(l) + "_" + to_string(Is3)}
    {
      // Retrieve CKM matrix element.
      const int i = Vij.at(l).first;
      const int j = Vij.at(l).second;

      _rules[GLUON]   = { {CG, GLUON, fact} };
      _rules[SIGMA]   = { {CS, SIGMA, fact / 6} };
      _rules[VALENCE] = { {CS, VALENCE, fact / 6} };
      for (int k = 2; k <= 6; k++)
	{
	  double coefp = 0;
	  double coefm = 0;
	  if (k >= 2 * j - 1)
	    {
	      double a = 1;
	      if (k == 2 * j - 1)
		a -= k;

	      coefp += a;
	      coefm += a;
	    }
	  if (k >= 2 * i)
	    {
	      double b = 1;
	      if (k == 2 * i)
		b -= k;

	      coefp += b;
	      coefm -= b;
	    }

	  coefp /= 2 * k * ( k - 1 );
	  coefm /= 2 * k * ( k - 1 );

	  // Change sign to T3 and V3
	  if (k == 2)
	    {
	      coefp *= - 1;
	      coefm *= - 1;
	    }

	  // Multiply the coefficient by the overall factor.
	  coefp *= fact;
	  coefm *= fact;

	  if (Is3)
	    {
	      _rules[2*k-1] = { {CNS, 2*k-1, coefm} };
	      _rules[2*k]   = { {CNS, 2*k,   coefp} };
	    }
	  else
	    {
	      _rules[2*k-1] = { {CNS, 2*k-1, coefp} };
	      _rules[2*k]   = { {CNS, 2*k,   coefm} };
	    }
	}
    };

    /**
     * @brief The DISCCBasis default constructor for the total structure functions.
     * @param CKM: vector with the CKM matrix entries
     * @param Is3: switch to tell the constructure whether the structure function is F3 or not
     */
  DISCCBasis(vector<double> const& CKM, bool const& Is3):
    ConvolutionMap{"DISCCBasis_tot"}
    {
      if (CKM.size() != 9)
	throw runtime_error(error("DISCCBasis", "The CKM vector must have 9 entries."));

      _rules[GLUON]   = { {CG, GLUON, 0} };
      _rules[SIGMA]   = { {CS, SIGMA, 0} };
      _rules[VALENCE] = { {CS, VALENCE, 0} };
      for (int k = 3; k <= 12; k++)
	_rules[k] = { {CNS, k, 0} };

      // Fill in rules according to the input CKM matrix.
      for (auto i = 1; i <= (int) CKM.size(); i++)
	{
	  if (CKM[i-1] == 0)
	    continue;

	  // Get basis of the i-th component.
	  const DISCCBasis basis{i, Is3};

	  // Get the rules.
	  const auto rules = basis.GetRules();

	  // Update rules
	  for (int k = 0; k <= 12; k++)
	    if (rules.count(k) != 0)
	      _rules[k][0].coefficient += CKM[i-1] * rules.at(k)[0].coefficient;
	}
    };
    ///@}
  };
  ///@}
}
