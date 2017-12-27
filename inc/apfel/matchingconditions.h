//
// APFEL++ 2017
//
// Authors: Valerio Bertone: valerio.bertone@cern.ch
//          Stefano Carrazza: stefano.carrazza@cern.ch
//

#pragma once

#include "apfel/expression.h"

namespace apfel
{
  /**
   * @notes Expressions taken from:
   * https://arxiv.org/pdf/hep-ph/9612398.pdf. Note that in these
   * expressions ln(m2/mu2) appears while we need ln(mu2/m2), so we
   * need to include a sign minus in front of every term linear in
   * this log.
   */

  /**
   * @brief O(alpha_s) matching conditions
   */
  // Term propotional to ln(mu2/m2) of eq. (B.2)
  class AS1Hg_L: public Expression
  {
  public:
    AS1Hg_L();
    double Regular(double const& x)  const;
  };

  // Term propotional to ln(mu2/m2) of eq (B.6)
  class AS1ggH_L: public Expression
  {
  public:
    AS1ggH_L();
    double Local(double const&)  const;
  };

  /**
   * @brief O(alpha_s^2) matching conditions
   */
  // Constant term of eq (B.1)
  class APS2Hq_0: public Expression
  {
  public:
    APS2Hq_0();
    double Regular(double const& x)  const;
  };

  // Constant term of eq (B.3)
  class AS2Hg_0: public Expression
  {
  public:
    AS2Hg_0();
    double Regular(double const& x)  const;
  };

  // Constant term of eq (B.4)
  class ANS2qqH_0: public Expression
  {
  public:
    ANS2qqH_0();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  // Constant term of eq (B.5)
  class AS2gqH_0: public Expression
  {
  public:
    AS2gqH_0();
    double Regular(double const& x)  const;
  };

  // Constant term of eq (B.7)
  class AS2ggH_0: public Expression
  {
  public:
    AS2ggH_0();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };
}
