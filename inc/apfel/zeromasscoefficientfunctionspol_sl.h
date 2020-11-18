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
   * @defgroup NCMasslesspol Zero-mass coefficient functions for longitudinally polarised DIS
   * Collection of the Zero-mass coefficient functions for
   * g<SUB>4</SUB>, g<SUB>L</SUB>, and g<SUB>1</SUB> up to
   * O(&alpha;<SUB>s</SUB>).
   */
  ///@{
  /**
   * @defgroup NLOzmpol NLO zero-mass coefficient functions
   * @ingroup NCMasslesspol
   */
  ///@{
  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for g<SUB>4</SUB>.
   */
  class G41ns: public Expression
  {
  public:
    G41ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for g<SUB>L</SUB>.
   */
  class GL1ns: public Expression
  {
  public:
    GL1ns();
    double Regular(double const& x) const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) non-singlet coefficient function
   * for g<SUB>1</SUB>.
   */
  class G11ns: public Expression
  {
  public:
    G11ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double Local(double const& x)    const;
  };

  /**
   * @brief O(&alpha;<SUB>s</SUB>) gluon coefficient function for for
   * g<SUB>1</SUB>.
   */
  class G11g: public Expression
  {
  public:
    G11g();
    double Regular(double const& x) const;
  };
  ///@}
  ///@}
}
