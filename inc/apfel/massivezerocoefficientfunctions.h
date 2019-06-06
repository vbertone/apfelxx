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
   * @defgroup NCMassiveZero Massless limit of the massive neutral current coefficient functions
   * Collection of the neutral current massless limit of the massive
   * coefficient functions for F<SUB>2</SUB> and F<SUB>L</SUB> up to
   * O(&alpha;<SUB>s</SUB>).
   * @note In the following "xi"
   * indicates the ratio Q<SUP>2</SUP> / M<SUP>2</SUP> and "xiF" is
   * equal to M<SUP>2</SUP> / &mu;<SUB>F</SUB><SUP>2</SUP>. See Appendix D of
   * https://arxiv.org/pdf/hep-ph/9601302.pdf. The suffix "_c" in the
   * classes below stands for "constant", that is no dependence on
   * "xi" or "xiF". "_l" or "_l2" instead means that that term is
   * proportional to log(xi) or log<SUP>2</SUP>(xi),
   * respectively. While "_f" means proportional to log(xiF) amd "_lf"
   * proportional to log(xi)log(xiF).
   */
  ///@{
  /**
   * @defgroup LOZero LO massive-zero coefficient functions
   * @ingroup NCMassiveZero
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for
   * F2. Constant term.
   */
  class Cm021gNC_c: public Expression
  {
  public:
    Cm021gNC_c();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for
   * F2. Single-log term.
   */
  class Cm021gNC_l: public Expression
  {
  public:
    Cm021gNC_l();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for
   * FL. Constant term.
   */
  class Cm0L1gNC_c: public Expression
  {
  public:
    Cm0L1gNC_c();
    double Regular(double const& x)  const;
  };
  ///@}

  /**
   * @defgroup NLOZero NLO massive-zero coefficient functions
   * @ingroup NCMassiveZero
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for F2. Constant term.
   */
  class Cm022nsNC_c: public Expression
  {
  public:
    Cm022nsNC_c();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for F2. Single-log term.
   */
  class Cm022nsNC_l: public Expression
  {
  public:
    Cm022nsNC_l();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for F2. Double-log term.
   */
  class Cm022nsNC_l2: public Expression
  {
  public:
    Cm022nsNC_l2();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for FL. Constant term.
   */
  class Cm0L2nsNC_c: public Expression
  {
  public:
    Cm0L2nsNC_c();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) non-singlet coefficient
   * function for FL. Single-log term.
   */
  class Cm0L2nsNC_l: public Expression
  {
  public:
    Cm0L2nsNC_l();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2. Constant term.
   */
  class Cm022psNC_c: public Expression
  {
  public:
    Cm022psNC_c();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2. Single-log term.
   */
  class Cm022psNC_l: public Expression
  {
  public:
    Cm022psNC_l();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2. Double-log term.
   */
  class Cm022psNC_l2: public Expression
  {
  public:
    Cm022psNC_l2();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2. Single-log(&mu;<SUB>F</SUB>) term.
   */
  class Cm022psNC_f: public Expression
  {
  public:
    Cm022psNC_f();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for F2. Mixed-double-log term.
   */
  class Cm022psNC_lf: public Expression
  {
  public:
    Cm022psNC_lf();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for FL. Constant term.
   */
  class Cm0L2psNC_c: public Expression
  {
  public:
    Cm0L2psNC_c();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for FL. Single-log term.
   */
  class Cm0L2psNC_l: public Expression
  {
  public:
    Cm0L2psNC_l();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) pure-singlet
   * coefficient function for FL. Single-log(&mu;<SUB>F</SUB>) term.
   */
  class Cm0L2psNC_f: public Expression
  {
  public:
    Cm0L2psNC_f();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2. Constant term.
   */
  class Cm022gNC_c: public Expression
  {
  public:
    Cm022gNC_c();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2. Single-log term.
   */
  class Cm022gNC_l: public Expression
  {
  public:
    Cm022gNC_l();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2. Double-log term.
   */
  class Cm022gNC_l2: public Expression
  {
  public:
    Cm022gNC_l2();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2. Single-log(&mu;<SUB>F</SUB>) term.
   */
  class Cm022gNC_f: public Expression
  {
  public:
    Cm022gNC_f();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for F2. Mixed-double-log term.
   */
  class Cm022gNC_lf: public Expression
  {
  public:
    Cm022gNC_lf();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for FL. Constant term.
   */
  class Cm0L2gNC_c: public Expression
  {
  public:
    Cm0L2gNC_c();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for FL. Single-log term.
   */
  class Cm0L2gNC_l: public Expression
  {
  public:
    Cm0L2gNC_l();
    double Regular(double const& x)  const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>2</SUP>) gluon coefficient
   * function for FL. Single-log(&mu;<SUB>F</SUB>) term.
   */
  class Cm0L2gNC_f: public Expression
  {
  public:
    Cm0L2gNC_f();
    double Regular(double const& x)  const;
  };
  ///@}
  ///@}
}
