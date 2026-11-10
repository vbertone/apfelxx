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
   * @defgroup CFFMassless Zero-mass coefficient functions for the Compton form factors
   * Collection of the MSbar coefficient functions for the compton
   * form factors.
   */
  ///@{
  ///@}

  /**
   * @defgroup CFFLOzm LO zero-mass coefficient functions
   * @ingroup CFFMassless
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB><SUP>0</SUP>) non-singlet real part of the
   * coefficient function for F1, identical for both signs of
   * Q<SUP>2</SUP>. The imaginary part is proportional to the identity
   * and the coefficient functions for FL vanish at this order.
   */
  class ReCFF10ns: public Expression
  {
  public:
    ReCFF10ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double LocalPP(double const& c, double const& d) const;
  };
  ///@}

  /**
   * @defgroup CFFNLOzm NLO zero-mass coefficient functions
   * @ingroup CFFMassless
   */
  ///@{

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet imaginary part of the
   * coefficient function for FL for Q<SUP>2</SUP> &gt; 0.
   */
  class ImCFFL1nsSpaceLike: public Expression
  {
  public:
    ImCFFL1nsSpaceLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet real part of the
   * coefficient function for FL for Q<SUP>2</SUP> &gt; 0.
   */
  class ReCFFL1nsSpaceLike: public Expression
  {
  public:
    ReCFFL1nsSpaceLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet imaginary part of the
   * coefficient function for FL for Q<SUP>2</SUP> &lt; 0.
   */
  class ImCFFL1nsTimeLike: public Expression
  {
  public:
    ImCFFL1nsTimeLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet real part of the
   * coefficient function for FL for Q<SUP>2</SUP> &lt; 0.
   */
  class ReCFFL1nsTimeLike: public Expression
  {
  public:
    ReCFFL1nsTimeLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon imaginary part of the
   * coefficient function for FL for Q<SUP>2</SUP> &gt; 0.
   */
  class ImCFFL1gSpaceLike: public Expression
  {
  public:
    ImCFFL1gSpaceLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon real part of the
   * coefficient function for FL for Q<SUP>2</SUP> &gt; 0.
   */
  class ReCFFL1gSpaceLike: public Expression
  {
  public:
    ReCFFL1gSpaceLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon imaginary part of the
   * coefficient function for FL for Q<SUP>2</SUP> &lt; 0.
   */
  class ImCFFL1gTimeLike: public Expression
  {
  public:
    ImCFFL1gTimeLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon real part of the
   * coefficient function for FL for Q<SUP>2</SUP> &lt; 0.
   */
  class ReCFFL1gTimeLike: public Expression
  {
  public:
    ReCFFL1gTimeLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet imaginary part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &gt; 0.
   */
  class ImCFF11nsSpaceLike: public Expression
  {
  public:
    ImCFF11nsSpaceLike(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& c)    const;
    double LocalPP(double const& c, double const& d) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet real part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &gt; 0.
   */
  class ReCFF11nsSpaceLike: public Expression
  {
  public:
    ReCFF11nsSpaceLike(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& c)    const;
    double LocalPP(double const& c, double const& d) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet imaginary part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &lt; 0.
   */
  class ImCFF11nsTimeLike: public Expression
  {
  public:
    ImCFF11nsTimeLike(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& c)    const;
    double LocalPP(double const& c, double const& d) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet real part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &lt; 0.
   */
  class ReCFF11nsTimeLike: public Expression
  {
  public:
    ReCFF11nsTimeLike(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& c)    const;
    double LocalPP(double const& c, double const& d) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon imaginary part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &gt; 0.
   */
  class ImCFF11gSpaceLike: public Expression
  {
  public:
    ImCFF11gSpaceLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon real part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &gt; 0.
   */
  class ReCFF11gSpaceLike: public Expression
  {
  public:
    ReCFF11gSpaceLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon imaginary part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &lt; 0.
   */
  class ImCFF11gTimeLike: public Expression
  {
  public:
    ImCFF11gTimeLike(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon real part of the
   * coefficient function for F1 for Q<SUP>2</SUP> &lt; 0.
   */
  class ReCFF11gTimeLike: public Expression
  {
  public:
    ReCFF11gTimeLike(double const& xi);
    double Regular(double const& y) const;
  };
  ///@}
}
