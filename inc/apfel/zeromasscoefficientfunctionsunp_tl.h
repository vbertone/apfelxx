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
   * @defgroup NCMasslessSIA Zero-mass coefficient functions for unpolarised SIA
   * Collection of the Zero-mass coefficient functions for
   * F<SUB>2</SUB>, F<SUB>L</SUB>, and F<SUB>3</SUB> up to
   * O(&alpha;<SUB>s</SUB><SUP>2</SUP>).
   */
  ///@{
  /**
   * @defgroup NLOzmSIA NLO zero-mass coefficient functions
   * @ingroup NCMasslessSIA
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for F2 in SIA.
   */
  class C21Tns: public Expression
  {
  public:
    C21Tns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for F2
   * in SIA.
   */
  class C21Tg: public Expression
  {
  public:
    C21Tg();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for FL in SIA.
   */
  class CL1Tns: public Expression
  {
  public:
    CL1Tns();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for FL
   * in SIA.
   */
  class CL1Tg: public Expression
  {
  public:
    CL1Tg();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for F3 in SIA.
   */
  class C31Tns: public Expression
  {
  public:
    C31Tns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
  ///@}

  /**
   * @defgroup NNLOzmSIA NNLO zero-mass coefficient functions
   * @ingroup NCMasslessSIA
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for F2 in SIA.
   */
  class C22Tnsp: public Expression
  {
  public:
    C22Tnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double _a0;
    double _a1;
    double _a2;
    double _a3;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2 in SIA.
   */
  class C22Tps: public Expression
  {
  public:
    C22Tps();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2 in SIA.
   */
  class C22Tg: public Expression
  {
  public:
    C22Tg();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for FL in SIA.
   */
  class CL2Tnsp: public Expression
  {
  public:
    CL2Tnsp(int const& nf);
    double Regular(double const& x) const;
  private:
    int const _nf;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for FL in SIA.
   */
  class CL2Tps: public Expression
  {
  public:
    CL2Tps();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for FL in SIA.
   */
  class CL2Tg: public Expression
  {
  public:
    CL2Tg();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet-plus
   * coefficient function for F3 in SIA.
   */
  class C32Tnsp: public Expression
  {
  public:
    C32Tnsp(int const& nf);
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  private:
    int const _nf;
    double _a0;
    double _a1;
    double _a2;
    double _a3;
  };
  ///@}
  ///@}
}
