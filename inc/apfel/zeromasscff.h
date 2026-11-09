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
   * coefficient function for F2, identical for both signs of
   * Q<SUP>2</SUP>. The imaginary part is proportional to the identity
   * and the coefficient functions for FL vanish at this order.
   */
  class ReCFF20ns: public Expression
  {
  public:
    ReCFF20ns();
    double Regular(double const& x)  const;
    double Singular(double const& x) const;
    double LocalPP(double const& c, double const& d) const;
  };
  ///@}
}
