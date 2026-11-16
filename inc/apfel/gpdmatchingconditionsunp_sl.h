//
// APFEL++ 2017
//
// Author: Valerio Bertone: valerio.bertone@cern.ch
//

#pragma once

#include "apfel/gpdsplittingfunctionsunp_sl.h"

namespace apfel
{
  /**
   * @defgroup GPDMatchCond GPD matching conditions
   * Collection of the MSbar matching conditions of unpolarised GPDs
   * at the heavy-quark thresholds. The terms that do not multiply
   * ln(&mu;<SUP>2</SUP>/m<SUB>h</SUB><SUP>2</SUP>) are only known
   * for parton-level skewness &xi;' = &kappa; y < 1 and are set to
   * zero for &xi;' > 1.
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term of the
   * gluon-in-heavy-quark matching function M<SUB>gh</SUB>. The
   * logarithmic term is given by Pgpd0gq.
   */
  class AS1gpdgH_0: public Expression
  {
  public:
    AS1gpdgH_0(double const& xi);
    double Regular(double const& y) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) constant term of the
   * heavy-quark-in-heavy-quark matching function M<SUB>hh</SUB>.
   */
  class AS1gpdHH_0: public Expression
  {
  public:
    AS1gpdHH_0(double const& xi);
    double Regular(double const& y)  const;
    double Singular(double const& y) const;
    double Local(double const& y)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) coefficient of
   * ln(&mu;<SUP>2</SUP>/m<SUB>h</SUB><SUP>2</SUP>) in the
   * heavy-quark-in-heavy-quark matching function M<SUB>hh</SUB>:
   * Pgpd0qq plus the cusp term 2 C<SUB>F</SUB>(&chi; coth &chi; - 1)
   * &delta;(1 - y).
   */
  class AS1gpdHH_L: public Pgpd0qq
  {
  public:
    AS1gpdHH_L(double const& xi);
    double Local(double const& y) const;
  };
  ///@}
}
