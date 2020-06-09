//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"

namespace apfel
{
  /**
   * @defgroup DISBases DIS convolution maps
   * Collection of derived classes from ConvolutionMap that implement
   * the convolution map for NC and CC DIS structure functions.
   * @ingroup ConvMap
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
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {CNS, CS, CG};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DISNCBasis constructor for the k-th
     * component of the structure functions.
     * @param k: index that identifies the component of the structure function
     * @param fact: factor that multiplies the whole structure function
     */
    DISNCBasis(int const& k, double const& fact = 1);

    /**
     * @brief The DISNCBasis constructor for the total structure functions.
     * @param Ch: vector of the effective quark charges
     */
    DISNCBasis(std::vector<double> const& Ch);
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
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {CNS, CS, CG};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

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
    std::map<int, std::pair<int, int>> Vij =
    {
      {1, {1, 1}}, {2, {1, 2}}, {3, {1, 3}},
      {4, {2, 1}}, {5, {2, 2}}, {6, {2, 3}},
      {7, {3, 1}}, {8, {3, 2}}, {9, {3, 3}}
    };

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DISCCBasis constructor for the (i,k)-th
     * component of the structure functions.
     * @param l: index that identifies the component of the structure function
     * @param Is3: switch to tell the constructure whether the structure function is F3 or not
     * @param fact: factor that multiplies the whole structure function
     */
    DISCCBasis(int const& l, bool const& Is3, double const& fact = 1);

    /**
     * @brief The DISCCBasis constructor for the total structure functions.
     * @param CKM: vector with the CKM matrix entries
     * @param Is3: switch to tell the constructure whether the structure function is F3 or not
     */
    DISCCBasis(std::vector<double> const& CKM, bool const& Is3);
    ///@}
  };
  ///@}
}
