#pragma once

#include "apfel/expression.h"

namespace apfel
{

  /**
   * @defgroup ACOT Massiv coefficient functions
   * Collection of massive coefficient functions for the
   * (full) ACOT scheme up to O((&alpha;<SUB>s</SUB>).
   * @note Only neutral current is available.
   */

  ///@{
  /**
   * @brief O(1) quark coefficient function for F2.
   * Eq. (2) of https://arxiv.org/pdf/hep-ph/9805233.
   */
  class Cm20qNC_ACOT: public Expression{
    public:
      Cm20qNC_ACOT(double const& eta);
      double Local(double const& x) const;

    private:
      double _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for F2.
   * Eq. (2.1) of https://doi.org/10.1007/BF01584394.
   */
  class Cm21gNC_ACOT: public Expression
  {
  public:
    Cm21gNC_ACOT(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon subtraction function for F2.
   * Eq. (16) of https://arxiv.org/pdf/hep-ph/9805233.
   */
  class Cm21gNC_sub_ACOT: public Expression
  {
  public:
    Cm21gNC_sub_ACOT(double const& eta);
    double Regular(double const& x) const;
  private:
    double _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for F2.
   * Eq. (8) of https://arxiv.org/pdf/hep-ph/9805233.
   */
  class Cm21qNC_ACOT: public Expression
  {
  public:
    Cm21qNC_ACOT(double const& eta);
    double Local(double const& x) const;
    double Singular(double const& x) const;
  private:
    double g(double const& x) const;
    double _sud;
    double _V2;
    double _S2;
    double _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark subtraction function for F2.
   * Eq. (12) of https://arxiv.org/pdf/hep-ph/9805233.
   */
  class Cm21qNC_sub_ACOT: public Expression
  {
  public:
    Cm21qNC_sub_ACOT(double const& eta);
    double Local(double const& x) const;
    double Singular(double const& x) const;
  private:
    double _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for FL.
   * Eq. (2.1) of https://doi.org/10.1007/BF01584394.
   */
  class CmL1gNC_ACOT : public Expression
  {
  public:
    CmL1gNC_ACOT(double const& eta);
    double Regular(double const& x)  const;
  };
  ///@}

  /**
   * @defgroup SACOT-chi Massiv coefficient functions
   * Collection of massive coefficient functions for the
   * SACOT-chi scheme up to O((&alpha;<SUB>s</SUB>).
   */
  ///@{

  /**
   * @brief O(1) quark coefficient function for F2.
   * Eq. (27) of https://arxiv.org/abs/hep-ph/0003035.
   */
  class Cm20qNC_ACOT_chi: public Expression{
    public:
      Cm20qNC_ACOT_chi(double const& eta);
      double Local(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon subtraction function for F2.
   * Eq. (27) of https://arxiv.org/abs/hep-ph/0003035.
   */
  class Cm21gNC_sub_ACOT_chi: public Expression
  {
  public:
    Cm21gNC_sub_ACOT_chi(double const& eta);
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for F2.
   * Eq. (2.1) of https://doi.org/10.1007/BF01584394.
   */
  class Cm21qNC_ACOT_chi : public Expression
  {
  public:
    Cm21qNC_ACOT_chi(double const& eta);
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for F3.
   * Eq. (8) of https://arxiv.org/pdf/hep-ph/9805233.
   */
  class Cm31qNC_ACOT_chi : public Expression
  {
  public:
    Cm31qNC_ACOT_chi(double const& eta);
    double Singular(double const& x) const;
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) quark coefficient function for FL.
   * Eq. (8) of https://arxiv.org/pdf/hep-ph/9805233.
   */
  class CmL1qNC_ACOT_chi : public Expression
  {
  public:
    CmL1qNC_ACOT_chi(double const& eta);
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon subtraction function for F2.
   * Eq. (27) of https://arxiv.org/abs/hep-ph/0003035.
   */
  class Cm21gCC_sub_ACOT: public Expression
  {
  public:
    Cm21gCC_sub_ACOT(double const& eta, double const& xi);
    double Regular(double const& x)  const;
  private:
    const double _xi;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for F2.
   * Eq. (2.1) of https://doi.org/10.1007/BF01584394.
   */
  class Cm21gCC_general_mass: public Expression
  {
  public:
    Cm21gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2);
    double Regular(double const& x)  const;
  private:
    const double _xi;
    const double _m1;
    const double _m2;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for FL.
   * Eq. (2.1) of https://doi.org/10.1007/BF01584394.
   */
  class CmL1gCC_general_mass: public Expression 
  {
  public:
    CmL1gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2);
    double Regular(double const& x)  const;
  private:
    const double _xi;
    const double _m1;
    const double _m2;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for F3.
   * Eq. (2.1) of https://doi.org/10.1007/BF01584394.
   */
  class Cm31gCC_general_mass: public Expression 
  {
  public:
    Cm31gCC_general_mass(double const& eta, double const& xi, double const& m1, double const& m2);
    double Regular(double const& x)  const;
  private:
    const double _xi;
    const double _m1;
    const double _m2;
  };
  ///@}

  /**
   * @defgroup aSACOT-chi Massiv coefficient functions
   * Collection of massive coefficient functions for the
   * approximative SACOT-chi scheme only for 
   * O((&alpha;<SUB>s</SUB><SUP>2</SUP>).
   * @note All of these are taken from 
   * https://www.liverpool.ac.uk/~avogt/coeff.html, where the 
   * classes ending in _0 are the parts of the coefficient, where
   * n<SUB>f</SUB>=0 and the classes ending in _nf are the parts
   * proportional to nf.
   */
  ///@{
  class C22nsm_aSACOT_chi_0 : public Expression
  {
  public:
    C22nsm_aSACOT_chi_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22nsm_aSACOT_chi_nf : public Expression
  {
  public:
    C22nsm_aSACOT_chi_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22nsp_aSACOT_chi_0 : public Expression
  {
  public:
    C22nsp_aSACOT_chi_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22nsp_aSACOT_chi_nf : public Expression
  {
  public:
    C22nsp_aSACOT_chi_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C22ps_aSACOT_chi: public Expression
  {
  public:
    C22ps_aSACOT_chi(double const& eta, bool const& xdependent);
    double Regular(double const& x) const;
  };

  class C22g_aSACOT_chi: public Expression
  {
  public:
    C22g_aSACOT_chi(double const& eta, bool const& xdependent);
    double Regular(double const& x) const;
    double Local(double const& x)   const;
  };

  class CL2nsm_aSACOT_chi_0 : public Expression
  {
  public:
    CL2nsm_aSACOT_chi_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  };

  class CL2nsm_aSACOT_chi_nf : public Expression
  {
  public:
    CL2nsm_aSACOT_chi_nf(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class CL2nsp_aSACOT_chi_0 : public Expression
  {
  public:
    CL2nsp_aSACOT_chi_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
  };

  class CL2nsp_aSACOT_chi_nf : public Expression
  {
  public:
    CL2nsp_aSACOT_chi_nf(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class CL2g_aSACOT_chi : public Expression
  {
  public:
    CL2g_aSACOT_chi(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class CL2ps_aSACOT_chi : public Expression
  {
  public:
    CL2ps_aSACOT_chi(double const& eta,bool const& xdependent);
    double Regular(double const& x)  const;
  };

  class C32nsm_aSACOT_chi_0 : public Expression
  {
  public:
    C32nsm_aSACOT_chi_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C32nsm_aSACOT_chi_nf : public Expression
  {
  public:
    C32nsm_aSACOT_chi_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C32nsp_aSACOT_chi_0 : public Expression
  {
  public:
    C32nsp_aSACOT_chi_0(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };

  class C32nsp_aSACOT_chi_nf : public Expression
  {
  public:
    C32nsp_aSACOT_chi_nf(double const& eta,bool const& xdependent);
    double Local(double const& x)    const;
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
  };
  ///@}
}