//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/convolutionmap.h"
#include "apfel/operator.h"

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

  /**
   * @brief The DISNCBasis_ACOT class is a derived of ConvolutionMap
   * specialised for the computation of the NC DIS structure
   * functions in the ACOT scheme.
   */
  class DISNCBasis_ACOT: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Operand: int {OGLUON, OSIGMA, OVALENCE, OT3, OV3, OT8, OV8, OT15, OV15, OT24, OV24, OT35, OV35};
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DISNCBasis_ACOT constructor for the total structure functions.
     * @param Ch: vector of the effective quark charges
     */
    DISNCBasis_ACOT(std::vector<double> const& Ch);
    ///@}

    /**
     * @brief Computes the change of basis from the physics basis to the QCD-evolution basis for light structure function.  
     * 
     * @param isPV switches between convolution basis for parity violating charges (F3) or not (F1,F2,FL)
     * @param gluon Gluon operators for each order
     * @param ns NS operators for each order AND each number of flavours -> you need only nf=3
     * @param ps PS operators for each order 
     * @return Returns the QCD evolution map in the order of C0,C1,C2
     */
    std::vector<std::map<int,Operator>> get_light_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps);
    /**
     * @brief Computes the change of basis from the physics basis to the QCD-evolution basis for charm structure function.  
     * 
     * @param isPV switches between convolution basis for parity violating charges (F3) or not (F1,F2,FL)
     * @param gluon Gluon operators for each order
     * @param ns NS operators for each order AND each number of flavours -> you need only nf=3,4
     * @param ps PS operators for each order 
     * @return Returns the QCD evolution map in the order of C0,C1,C2
     */
    std::vector<std::map<int,Operator>> get_charm_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps);
    /**
     * @brief Computes the change of basis from the physics basis to the QCD-evolution basis for bottom structure function.  
     * 
     * @param isPV switches between convolution basis for parity violating charges (F3) or not (F1,F2,FL)
     * @param gluon Gluon operators for each order
     * @param ns NS operators for each order AND each number of flavours -> you need only nf=4,5
     * @param ps PS operators for each order 
     * @return Returns the QCD evolution map in the order of C0,C1,C2
     */
    std::vector<std::map<int,Operator>> get_bottom_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps);
    /**
     * @brief Computes the change of basis from the physics basis to the QCD-evolution basis for top structure function.  
     * 
     * @param isPV switches between convolution basis for parity violating charges (F3) or not (F1,F2,FL)
     * @param gluon Gluon operators for each order
     * @param ns NS operators for each order AND each number of flavours -> you need only nf=5,6
     * @param ps PS operators for each order 
     * @return Returns the QCD evolution map in the order of C0,C1,C2
     */
    std::vector<std::map<int,Operator>> get_top_operators(bool isPV, std::vector<Operator> gluon, std::vector<std::map<int,Operator>> ns, std::vector<Operator> ps);
    /**
     * @brief Uses the results from the light,charm,bottom and top get_operators-functions to calculate the total structure function. (It simply adds all coefficients for all orders and distributions)
     * 
     * @param isPV switches between convolution basis for parity violating charges (F3) or not (F1,F2,FL)
     * @param coeff The results from the light,charm,bottom and top get_operators-funcitons. Order: {light,charm,bottom,top}
     * @return Returns the QCD evolution map in the order of C0,C1,C2
     */
    std::vector<std::map<int,Operator>> get_tot_operators(bool isPV, std::vector<std::vector<std::map<int,Operator>>> coeff);
  private:
    /**
     * @brief Vector of effective charges and vector of averages as a function of nf
     */
    std::vector<double> _Ch;
    std::vector<double> _avCh;
  };

  /**
   * @brief The DISCCBasis_ACOT class is a derived of ConvolutionMap
   * specialised for the computation of the CC DIS structure
   * functions in the ACOT scheme.
   */
  class DISCCBasis_ACOT: public ConvolutionMap
  {
  public:
    /**
     * @brief The map enumerators for the operands and the
     * distributions.
     */
    enum Object:  int {GLUON, SIGMA, VALENCE, T3, V3, T8, V8, T15, V15, T24, V24, T35, V35};
    enum PhysDist: int {DOWN,UP,STRANGE,CHARM,BOTTOM,TOP};

    /**
     * @name Constructors
     * List of constructors.
     */
    ///@{
    /**
     * @brief The DISCCBasis_ACOT constructor
     * @param CKM: vector with the squared CKM matrix entries.
     * @param Zero: zero operator for initialization purposes
     */
    DISCCBasis_ACOT(std::vector<double> const& CKM, Operator Zero);
    ///@}

    /**
     * @brief Computes the change of basis from the QCD evolution to physical basis for the plus-operators
     * 
     * @param op_map: The operator map in the physical basis
     */
    std::vector<std::map<int,Operator>> get_operators_plus(std::vector<std::map<int,Operator>> op_map);
    /**
     * @brief Computes the change of basis from the QCD evolution to physical basis for the minus-operators
     * 
     * @param op_map: The operator map in the physical basis
     */
    std::vector<std::map<int,Operator>> get_operators_minus(std::vector<std::map<int,Operator>> op_map);

  private:
    /**
     * @brief The CKM matrix and the Zero-Operator
     * 
     */
    std::vector<double> const& _CKM;
    Operator _Zero;
  };
  ///@}
}
