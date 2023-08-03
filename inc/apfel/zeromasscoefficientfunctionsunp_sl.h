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

  /**
   * @defgroup NNNLOzm NNNLO zero-mass coefficient functions
   * @ingroup NCMassless
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-plus
   * coefficient function for F2.
   */
  class C23nsp: public Expression
  {
  public:
    C23nsp(int const& nf, bool const& fl11 = true);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-minus
   * coefficient function for F2.
   */
  class C23nsm: public Expression
  {
  public:
    C23nsm(int const& nf, bool const& fl11 = true);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * coefficient function for F2.
   */
  class C23ps: public Expression
  {
  public:
    C23ps(int const& nf, bool const& fl11 = true);
    double Regular(double const& x) const;
    double Local(double const& x)   const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for F2.
   */
  class C23g: public Expression
  {
  public:
    C23g(int const& nf, bool const& fl11 = true);
    double Regular(double const& x) const;
    double Local(double const& x)   const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-plus
   * coefficient function for FL.
   */
  class CL3nsp: public Expression
  {
  public:
    CL3nsp(int const& nf, bool const& fl11 = true);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-minus
   * coefficient function for FL.
   */
  class CL3nsm: public Expression
  {
  public:
    CL3nsm(int const& nf, bool const& fl11 = true);
    double Regular(double const& x) const;
    double Local(double const&)     const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) pure-singlet
   * coefficient function for FL.
   */
  class CL3ps: public Expression
  {
  public:
    CL3ps(int const& nf, bool const& fl11 = true);
    double Regular(double const& x) const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) gluon coefficient
   * function for FL.
   */
  class CL3g: public Expression
  {
  public:
    CL3g(int const& nf, bool const& fl11 = true);
    double Regular(double const& x) const;
  private:
    int  const _nf;
    bool const _fl11;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-plus
   * coefficient function for F3.
   */
  class C33nsp: public Expression
  {
  public:
    C33nsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) non-singlet-minus
   * coefficient function for F3.
   */
  class C33nsm: public Expression
  {
  public:
    C33nsm(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>3</SUP>) total-valence
   * coefficient function for F3.
   */
  class C33nsv: public Expression
  {
  public:
    C33nsv();
    double Regular(double const& x) const;
  };
  ///@}
  ///@}
}
