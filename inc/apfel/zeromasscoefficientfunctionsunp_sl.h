//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @defgroup NCMassless Zero-mass coefficient functions for unpolarised DIS
   * Collection of the Zero-mass coefficient functions for
   * F<SUB>2</SUB>, F<SUB>L</SUB>, and F<SUB>3</SUB> up to
   * O(&alpha;<SUB>s</SUB><SUP>2</SUP>).
   * @note While for the O(&alpha;<SUB>s</SUB>) coefficient functions
   * exact expressions are used, a fast parameterisation for the
   * O(&alpha;<SUB>s</SUB><SUP>2</SUP>) ones is used. See
   * https://www.liverpool.ac.uk/~avogt/coeff.html for more details.
   */
  ///@{
  /**
   * @defgroup NLOzm NLO zero-mass coefficient functions
   * @ingroup NCMassless
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for F2.
   */
  class C21ns: public Expression
  {
  public:
    C21ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for F2.
   */
  class C21g: public Expression
  {
  public:
    C21g();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for FL.
   */
  class CL1ns: public Expression
  {
  public:
    CL1ns();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for FL.
   */
  class CL1g: public Expression
  {
  public:
    CL1g();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for F3.
   */
  class C31ns: public Expression
  {
  public:
    C31ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NNLOzm NNLO zero-mass coefficient functions
   * @ingroup NCMassless
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for F2.
   */
  class C22nsp: public Expression
  {
  public:
    C22nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-minus
   * coefficient function for F2.
   */
  class C22nsm: public Expression
  {
  public:
    C22nsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2.
   */
  class C22ps: public Expression
  {
  public:
    C22ps();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2.
   */
  class C22g: public Expression
  {
  public:
    C22g();
    double Regular(double const& x) const;
    double Local(double const& x)   const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for FL.
   */
  class CL2nsp: public Expression
  {
  public:
    CL2nsp(int const& nf);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-minus
   * coefficient function for FL.
   */
  class CL2nsm: public Expression
  {
  public:
    CL2nsm(int const& nf);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for FL.
   */
  class CL2ps: public Expression
  {
  public:
    CL2ps();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for FL.
   */
  class CL2g: public Expression
  {
  public:
    CL2g();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for F3.
   */
  class C32nsp: public Expression
  {
  public:
    C32nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-minus
   * coefficient function for F3.
   */
  class C32nsm: public Expression
  {
  public:
    C32nsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };
  ///@}
  ///@}
}
